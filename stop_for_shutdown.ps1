param([ValidateSet('reference','group')][string]$Case='reference',[ValidateSet('base','time')][string]$Variant='base',[int]$TargetCompleted=9)
$ErrorActionPreference='Stop'
$workspaceRoot=$PSScriptRoot
$env:TEMP='E:/comsol_runtime/transient_matched/queue_temp'
$env:TMP=$env:TEMP
$root=Join-Path $workspaceRoot 'output/transient_converged'
$verified=Join-Path $workspaceRoot 'output/transient_verified'
$statusFile=Join-Path $verified 'shutdown_stop_status.json'
$progressFile=Join-Path $root "bulk/$Variant/L8/seed_101/$Case/blocks/P/progress.json"
function Read-Json($p){[IO.File]::ReadAllText($p,[Text.Encoding]::UTF8) | ConvertFrom-Json}
function State($v){$v | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statusFile -Encoding UTF8}
function Owned($path,$script){
    $m=Read-Json $path
    $p=Get-CimInstance Win32_Process -Filter "ProcessId=$($m.launcher_pid)"
    if($null -eq $p -or $p.Name -ne 'python.exe' -or $p.CommandLine -notlike "*$script*"){throw "Process identity changed: $script"}
    if([math]::Abs((([datetimeoffset]$p.CreationDate)-[datetimeoffset]::Parse($m.started)).TotalSeconds) -gt 2){throw 'PID was reused'}
    @(Get-CimInstance Win32_Process -Filter "ParentProcessId=$($p.ProcessId)" | Where-Object {$_.Name -eq 'python.exe' -and $_.CommandLine -like "*$script*"})+@($p)
}
$mutex=New-Object System.Threading.Mutex($false,'Local\CFRPShutdownCheckpointStop')
if(-not $mutex.WaitOne(0)){Write-Output 'Shutdown checkpoint stop already armed';exit 0}
try {
    $primary=@(Owned (Join-Path $root 'pipeline_process.json') 'run_refined_pipeline.py')
    State @{stage='waiting_for_checkpoint';case=$Case;variant=$Variant;target_completed=$TargetCompleted;primary_pids=@($primary.ProcessId)}
    while((Read-Json $progressFile).completed -lt $TargetCompleted){
        $current=@(Owned (Join-Path $root 'pipeline_process.json') 'run_refined_pipeline.py')
        if(($current.ProcessId -join ',') -ne ($primary.ProcessId -join ',')){throw 'Primary changed; refusing unrelated process stop'}
        Start-Sleep -Seconds 2
    }
    if((Read-Json (Join-Path $verified 'continuation_state.json')).stage -ne 'waiting_for_bulk_transient_controls'){throw 'Continuation is not idle'}
    $dependent=@(Owned (Join-Path $verified 'continuation_process.json') 'continue_transient_validation.py')
    $watcher=@(Owned (Join-Path $root 'watcher_process.json') 'watch_refined_results.py')
    $all=@($dependent)+@($watcher)+@($primary)
    $guidedRoot=Join-Path $workspaceRoot 'output/detectability_200kHz_L48'
    $guidedMetadata=Join-Path $guidedRoot 'worker_process.json'
    if(Test-Path -LiteralPath $guidedMetadata){
        $guidedMeta=Read-Json $guidedMetadata
        if(Get-CimInstance Win32_Process -Filter "ProcessId=$($guidedMeta.launcher_pid)"){
            $guidedState=Read-Json (Join-Path $guidedRoot 'worker_state.json')
            if($guidedState.stage -ne 'waiting_resources'){throw '48 mm guided worker is not at a safe wait boundary'}
            $guided=@(Owned $guidedMetadata 'detectability_200kHz_L48.py')
            $all+=@($guided)
        }
    }
    State @{stage='stopping_after_committed_checkpoint';case=$Case;variant=$Variant;target_completed=$TargetCompleted;stopping_pids=@($all.ProcessId)}
    foreach($p in $all){Stop-Process -Id $p.ProcessId -ErrorAction SilentlyContinue}
    Start-Sleep -Seconds 2
    foreach($p in $all){if(Get-Process -Id $p.ProcessId -ErrorAction SilentlyContinue){throw 'Worker has not exited'}}
    & (Join-Path $workspaceRoot '.venv/Scripts/python.exe') (Join-Path $workspaceRoot 'verify_shutdown_checkpoint.py') --case $Case --variant $Variant --minimum-completed $TargetCompleted
    if($LASTEXITCODE -ne 0){throw 'Checkpoint verification failed; not ready for shutdown'}
    State @{stage='ready_for_shutdown';case=$Case;variant=$Variant;stopped_pids=@($all.ProcessId);verified=(Get-Date).ToString('o');handoff='SHUTDOWN_HANDOFF.md'}
    Write-Output 'READY FOR SHUTDOWN'
} catch {State @{stage='needs_attention';error=$_.Exception.Message;at=(Get-Date).ToString('o')};throw}
finally {$mutex.ReleaseMutex();$mutex.Dispose()}
