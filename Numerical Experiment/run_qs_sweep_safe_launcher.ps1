$ErrorActionPreference = 'Stop'

$matlab = 'D:\matlab\bin\matlab.exe'
$workDir = 'D:\Mine\Mani\ManiCSDP\Numerical Experiment'
$outDir = Join-Path $workDir 'results_qs_sweep_safe'
$matlabLog = Join-Path $outDir 'matlab_full.log'
$launcherLog = Join-Path $outDir 'launcher.log'

New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$batch = "addpath(genpath('D:\Mine\Mani\ManiCSDP')); addpath(genpath('D:\Mine\Mani\ManiSDP-matlab-main\ManiSDP-matlab-main')); addpath(genpath('D:\matlab\cvx\sedumi')); cd('$workDir'); run_qs_sweep_safe"

"Launcher started: $(Get-Date -Format o)" | Out-File -FilePath $launcherLog -Encoding utf8
"MATLAB: $matlab" | Out-File -FilePath $launcherLog -Encoding utf8 -Append
"MATLAB log: $matlabLog" | Out-File -FilePath $launcherLog -Encoding utf8 -Append
"Batch: $batch" | Out-File -FilePath $launcherLog -Encoding utf8 -Append

& $matlab -wait -logfile $matlabLog -batch $batch 2>&1 |
    Out-File -FilePath $launcherLog -Encoding utf8 -Append

"Launcher finished: $(Get-Date -Format o)" | Out-File -FilePath $launcherLog -Encoding utf8 -Append
"Exit code: $LASTEXITCODE" | Out-File -FilePath $launcherLog -Encoding utf8 -Append
