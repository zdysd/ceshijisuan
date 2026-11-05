param(
    [ValidateSet('science','comsol','study','transient','detectability','tests','doctor')][string]$Task = 'doctor',
    [Parameter(ValueFromRemainingArguments = $true)][string[]]$TaskArguments
)
$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    $porePython = if ($env:CFRP_COMSOL_PYTHON) { $env:CFRP_COMSOL_PYTHON } else { Join-Path $PSScriptRoot '.venv/Scripts/python.exe' }
    if (-not (Test-Path -LiteralPath $porePython)) { throw 'Run setup.ps1 or set CFRP_COMSOL_PYTHON to a Python 3.12 environment with the project dependencies.' }
    & $porePython tools/restore_inputs.py
    if ($LASTEXITCODE -ne 0) { throw 'Bundled input verification failed.' }
    switch ($Task) {
        'tests' { & $porePython -B -m unittest discover -s tests -v }
        'doctor' { & $porePython -B tools/check_environment.py @TaskArguments }
        'science' { & $porePython experiment.py @TaskArguments }
        'detectability' {
            if (-not $TaskArguments) { $TaskArguments = @('report') }
            & $porePython -u detectability_200kHz_L48.py @TaskArguments
        }
        default {
            $poreScript = switch ($Task) { 'comsol' {'run_comsol.py'} 'study' {'run_study.py'} 'transient' {'run_transient.py'} }
            & $porePython -u $poreScript @TaskArguments
        }
    }
    if ($LASTEXITCODE -ne 0) { throw "Task failed with exit code $LASTEXITCODE" }
} finally { Pop-Location }
