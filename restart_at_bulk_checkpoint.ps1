# One-time, guarded restart at the next committed bulk block. Never rerun the pulse.
param([switch]$StatusOnly)
$ErrorActionPreference='Stop'
$workspaceRoot=$PSScriptRoot
$root=Join-Path $workspaceRoot 'output/transient_converged'
$verified=Join-Path $workspaceRoot 'output/transient_verified'
$progressFile=Join-Path $root 'bulk/base/L8/seed_101/reference/blocks/P/progress.json'
$statusFile=Join-Path $root 'checkpoint_restart.json'
$mutex=New-Object System.Threading.Mutex($false,'Local\CFRPBulkCheckpointRestart')
if(-not $mutex.WaitOne(0)){Write-Output 'A checkpoint restart is already armed.';exit 0}
function Read-Json($path){Get-Content -LiteralPath $path -Raw | ConvertFrom-Json}
function Write-State($state){$state | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statusFile -Encoding UTF8}
function Owned-Processes($metadataPath,$scriptName){
    $metadata=Read-Json $metadataPath
    $launcher=Get-CimInstance Win32_Process -Filter "ProcessId=$($metadata.launcher_pid)"
    if($null -eq $launcher -or $launcher.Name -ne 'python.exe' -or $launcher.CommandLine -notlike "*$scriptName*"){throw 'Registered process identity changed'}
    $expected=[datetimeoffset]::Parse($metadata.started)
    if([math]::Abs((([datetimeoffset]$launcher.CreationDate)-$expected).TotalSeconds) -gt 2){throw 'Registered process ID was reused'}
    $children=@(Get-CimInstance Win32_Process -Filter "ParentProcessId=$($launcher.ProcessId)" | Where-Object {$_.Name -eq 'python.exe' -and $_.CommandLine -like "*$scriptName*"})
    @($children)+@($launcher)
}
try {
    $primary=@(Owned-Processes (Join-Path $root 'pipeline_process.json') 'run_refined_pipeline.py')
    $target=(Read-Json $progressFile).completed+1
    if($StatusOnly){@{registered_pids=@($primary.ProcessId);next_completed_block=$target} | ConvertTo-Json;return}
    Write-State @{stage='waiting_for_committed_block';target_completed=$target;registered_pids=@($primary.ProcessId);purpose='Load reference 50 us / group 100 us caps; keep bounded automatic doubling and committed histories'}
    do {
        Start-Sleep -Seconds 5
        $current=@(Owned-Processes (Join-Path $root 'pipeline_process.json') 'run_refined_pipeline.py')
        if(($current.ProcessId -join ',') -ne ($primary.ProcessId -join ',')){throw 'Primary process set changed; no restart performed'}
        $progress=Read-Json $progressFile
    } while($progress.completed -lt $target)
    $continuation=Read-Json (Join-Path $verified 'continuation_state.json')
    if($continuation.stage -ne 'waiting_for_bulk_transient_controls'){throw 'Continuation worker is not idle; leave all solvers running'}
    $dependent=@(Owned-Processes (Join-Path $verified 'continuation_process.json') 'continue_transient_validation.py')
    Write-State @{stage='restarting_after_commit';completed=$progress.completed;checkpoint_file=$progress.checkpoint_file;target_completed=$target}
    foreach($taskProcess in $dependent){Stop-Process -Id $taskProcess.ProcessId -ErrorAction SilentlyContinue}
    foreach($taskProcess in $primary){Stop-Process -Id $taskProcess.ProcessId -ErrorAction SilentlyContinue}
    Start-Sleep -Seconds 2
    & (Join-Path $workspaceRoot 'resume_transient.ps1')
    Write-State @{stage='resumed';completed=$progress.completed;target_completed=$target;finished=(Get-Date).ToString('o')}
} catch {
    Write-State @{stage='restart_not_completed';error=$_.Exception.Message;recorded=(Get-Date).ToString('o')}
    throw
} finally {$mutex.ReleaseMutex();$mutex.Dispose()}
