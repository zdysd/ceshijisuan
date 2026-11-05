param([switch]$StatusOnly,[switch]$RetryAfterFix,[switch]$ResumeAfterShutdown)
$ErrorActionPreference='Stop'
$workspaceRoot=$PSScriptRoot
$outputRoot=Join-Path $workspaceRoot 'output/transient_converged'
$verifiedRoot=Join-Path $workspaceRoot 'output/transient_verified'
$pythonLauncher='E:/MCP/venvs/comsol/Scripts/python.exe'
$runtimeRoot='E:/comsol_runtime/transient_matched/queue_temp'
$env:TEMP=$runtimeRoot
$env:TMP=$runtimeRoot
$env:MPLCONFIGDIR='E:/comsol_runtime/transient_matched/matplotlib'
New-Item -ItemType Directory -Path $runtimeRoot,$verifiedRoot -Force | Out-Null
$mutex=New-Object System.Threading.Mutex($false,'Local\CFRPTransientResume')
if(-not $mutex.WaitOne(0)){Write-Output 'Another resume command is already running.';exit 0}
try {
    $stopStatePath=Join-Path $verifiedRoot 'shutdown_stop_status.json'
    if(Test-Path -LiteralPath $stopStatePath){
        $stopState=[IO.File]::ReadAllText($stopStatePath,[Text.Encoding]::UTF8) | ConvertFrom-Json
        if($stopState.stage -in @('waiting_for_checkpoint','stopping_after_committed_checkpoint')){
            Write-Output 'User-requested checkpoint stop is pending. Leave the current workers and stop guard running.'
            return
        }
    }
    $pauseFile=Join-Path $verifiedRoot 'PAUSED_FOR_SHUTDOWN.json'
    if(Test-Path -LiteralPath $pauseFile){
        if(-not $ResumeAfterShutdown){Write-Output 'Paused for shutdown. Explicitly resume with -ResumeAfterShutdown.';return}
        if(-not $StatusOnly){Move-Item -LiteralPath $pauseFile -Destination (Join-Path $verifiedRoot ('shutdown_pause_'+(Get-Date -Format 'yyyyMMdd_HHmmss')+'.json'))}
    }
    function Read-Json($path){if(Test-Path -LiteralPath $path){[IO.File]::ReadAllText($path,[Text.Encoding]::UTF8) | ConvertFrom-Json}}
    function Script-Hash($path){
        $algorithm=[Security.Cryptography.SHA256]::Create()
        $stream=[IO.File]::OpenRead($path)
        try {[BitConverter]::ToString($algorithm.ComputeHash($stream)).Replace('-','')}
        finally {$stream.Dispose();$algorithm.Dispose()}
    }
    function Running($scriptName){
        @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" | Where-Object { $_.CommandLine -like "*$scriptName*" })
    }
    if(@(Running 'salvage_bulk_recovery.py').Count){
        Write-Output 'Stopped-recovery validation is active; wait for its checkpoint publication before resuming the solver.'
        exit 0
    }
    function Start-JobProcess($scriptName,$arguments,$logDirectory,$logPrefix,$metadataPath){
        if(@(Running $scriptName).Count){Write-Output "Already running: $scriptName";return}
        if($StatusOnly){Write-Output "Would resume: $scriptName";return}
        $scriptHash=Script-Hash (Join-Path $workspaceRoot $scriptName)
        $stamp=Get-Date -Format 'yyyyMMdd_HHmmss'
        foreach($suffix in @('stdout','stderr')){
            $path=Join-Path $logDirectory ($logPrefix+'_'+$suffix+'.log')
            if(Test-Path -LiteralPath $path){
                $resolved=[IO.Path]::GetFullPath($path)
                $destination=[IO.Path]::GetFullPath((Join-Path $logDirectory ($logPrefix+'_'+$suffix+'_'+$stamp+'.log')))
                $allowed=[IO.Path]::GetFullPath($workspaceRoot).TrimEnd('\')+'\'
                if(-not $resolved.StartsWith($allowed,[StringComparison]::OrdinalIgnoreCase) -or
                   -not $destination.StartsWith($allowed,[StringComparison]::OrdinalIgnoreCase)){throw 'Log archive is outside the workspace'}
                Move-Item -LiteralPath $resolved -Destination $destination
            }
        }
        $job=Start-Process -FilePath $pythonLauncher -ArgumentList (@('-u',$scriptName)+$arguments) -WorkingDirectory $workspaceRoot -WindowStyle Hidden -PassThru `
            -RedirectStandardOutput (Join-Path $logDirectory ($logPrefix+'_stdout.log')) `
            -RedirectStandardError (Join-Path $logDirectory ($logPrefix+'_stderr.log'))
        @{launcher_pid=$job.Id;started=$job.StartTime.ToString('o');script=$scriptName;arguments=$arguments;
          temporary_directory=$runtimeRoot;script_sha256=$scriptHash} |
          ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $metadataPath -Encoding UTF8
        Write-Output "Resumed $scriptName, launcher PID $($job.Id)"
    }
    $primaryNeeded=$false
    foreach($kind in @('guided','bulk')){foreach($variant in @('base','time')){foreach($case in @('reference','group')){
        $record=Read-Json (Join-Path $outputRoot "$kind/$variant/L8/seed_101/$case/run.json")
        if($null -eq $record -or $record.status -ne 'complete'){$primaryNeeded=$true}
    }}}
    if($primaryNeeded){Start-JobProcess 'run_refined_pipeline.py' @('all','--profile','refined_v6','--cores','4','--samples','12') $outputRoot 'pipeline' (Join-Path $outputRoot 'pipeline_process.json')}
    else {Write-Output 'Original transient cases are complete; their completion does not imply acceptance.'}
    $probe=Join-Path $workspaceRoot 'output/transient_refined/guided_sampling'
    if(-not (Test-Path -LiteralPath (Join-Path $probe 'sampling_validation.json'))){
        Start-JobProcess 'check_guided_sampling.py' @('--cores','2') $outputRoot 'sampling' (Join-Path $outputRoot 'sampling_process.json')
    }
    $continuation=Read-Json (Join-Path $verifiedRoot 'continuation_state.json')
    if($continuation.stage -eq 'needs_adjustment' -and -not $RetryAfterFix){
        Write-Output 'A scientific check failed. Inspect transient_verified/validation.json and continuation_state.json; use -RetryAfterFix only after addressing it.'
    } elseif($continuation.stage -ne 'complete'){
        Start-JobProcess 'continue_transient_validation.py' @() $verifiedRoot 'continuation' (Join-Path $verifiedRoot 'continuation_process.json')
    }
    Start-JobProcess 'watch_refined_results.py' @() $outputRoot 'watcher' (Join-Path $outputRoot 'watcher_process.json')
    if($ResumeAfterShutdown -and -not $StatusOnly){
        @{stage='resumed';resumed_at=(Get-Date).ToString('o');reason='User explicitly resumed after shutdown'} | ConvertTo-Json | Set-Content -LiteralPath $stopStatePath -Encoding UTF8
    }
    if($ResumeAfterShutdown){
        $guidedRoot=Join-Path $workspaceRoot 'output/detectability_200kHz_L48'
        $guidedState=Read-Json (Join-Path $guidedRoot 'worker_state.json')
        if($guidedState.stage -eq 'waiting_resources'){
            Start-JobProcess 'detectability_200kHz_L48.py' @('run') $guidedRoot 'worker' (Join-Path $guidedRoot 'worker_process.json')
        }
    }
    Write-Output 'Temporary and recovery files remain on E. Numerical failures never authorize statistical expansion.'
} finally {$mutex.ReleaseMutex();$mutex.Dispose()}
