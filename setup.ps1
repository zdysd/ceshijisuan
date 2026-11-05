param([string]$Python = 'py')
$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    if ($PSScriptRoot -match '[^\x00-\x7F]') { throw 'Clone this COMSOL/JPype project into an ASCII path, e.g. D:/Projects/cfrp.' }
    $porePython = Join-Path $PSScriptRoot '.venv/Scripts/python.exe'
    if (-not (Test-Path -LiteralPath $porePython)) {
        if ($Python -eq 'py') { & $Python -3.12 -m venv .venv } else { & $Python -m venv .venv }
        if ($LASTEXITCODE -ne 0) { throw 'Creating Python 3.12 environment failed.' }
    }
    & $porePython -c "import sys; assert sys.version_info[:2] == (3,12), 'Python 3.12 required'"
    if ($LASTEXITCODE -ne 0) { throw 'Python version mismatch.' }
    & $porePython -m pip install -r requirements-lock.txt
    if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
    & $porePython tools/restore_inputs.py
    if ($LASTEXITCODE -ne 0) { throw 'Input restoration failed.' }
    & $porePython -B tools/check_environment.py
    if ($LASTEXITCODE -ne 0) { throw 'Offline preflight failed.' }
    & $porePython -B -m unittest discover -s tests -v
    if ($LASTEXITCODE -ne 0) { throw 'Offline tests failed.' }
    Write-Output 'Offline setup passed. No COMSOL solve or license test was started.'
} finally { Pop-Location }
