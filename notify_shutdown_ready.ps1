# User-requested local notification, only after checkpoint verification and worker exit.
$ErrorActionPreference='Stop'
$root=$PSScriptRoot
$statusPath=Join-Path $root 'output/transient_verified/shutdown_stop_status.json'
$notificationPath=Join-Path $root 'output/transient_verified/shutdown_notification.json'
Add-Type -AssemblyName System.Windows.Forms
while($true){
    if(Test-Path -LiteralPath $statusPath){
        try {$state=Get-Content -LiteralPath $statusPath -Raw | ConvertFrom-Json}
        catch {Start-Sleep -Seconds 2;continue}
        if($state.stage -eq 'ready_for_shutdown'){
            $proof=Get-Content -LiteralPath (Join-Path $root 'output/transient_verified/PAUSED_FOR_SHUTDOWN.json') -Raw | ConvertFrom-Json
            if($proof.stage -ne 'ready_for_shutdown'){throw 'Shutdown verification marker missing'}
            $alive=@(Get-CimInstance Win32_Process -Filter "Name='python.exe'" | Where-Object {$_.CommandLine -match 'run_refined_pipeline.py|continue_transient_validation.py|watch_refined_results.py|detectability_200kHz_L48.py run'})
            if($alive.Count){throw 'Simulation worker is still running; do not show ready notification'}
            @{requested_at=(Get-Date).ToString('o');end_time_us=$proof.end_time_us;message='Verified checkpoint saved and workers stopped. Safe for normal shutdown.'} | ConvertTo-Json | Set-Content -LiteralPath $notificationPath -Encoding UTF8
            # Unicode escapes are unnecessary: load the UTF-8 handoff written by Python.
            $title=-join @([char]0x53ef,[char]0x4ee5,[char]0x5173,[char]0x673a)
            $body=Get-Content -LiteralPath (Join-Path $root 'output/transient_verified/shutdown_notification.txt') -Raw -Encoding UTF8
            [System.Windows.Forms.MessageBox]::Show($body,$title,[System.Windows.Forms.MessageBoxButtons]::OK,[System.Windows.Forms.MessageBoxIcon]::Information,[System.Windows.Forms.MessageBoxDefaultButton]::Button1,[System.Windows.Forms.MessageBoxOptions]::ServiceNotification) | Out-Null
            break
        }
        if($state.stage -eq 'needs_attention'){
            [System.Windows.Forms.MessageBox]::Show("Checkpoint stop needs attention. Do not assume shutdown preparation succeeded. See: $statusPath",'Simulation checkpoint check',[System.Windows.Forms.MessageBoxButtons]::OK,[System.Windows.Forms.MessageBoxIcon]::Warning) | Out-Null
            break
        }
    }
    Start-Sleep -Seconds 5
}
