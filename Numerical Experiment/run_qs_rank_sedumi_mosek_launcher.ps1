$ErrorActionPreference = "Continue"

$workDir = "D:\Mine\Mani\ManiCSDP\Numerical Experiment"
$outDir = Join-Path $workDir "results_qs_rank_sedumi_mosek"
$matlab = "D:\matlab\bin\matlab.exe"
$summary = Join-Path $outDir "launcher_summary.log"

New-Item -ItemType Directory -Force -Path $outDir | Out-Null

"QS rank sweep launcher started at $(Get-Date -Format o)" | Set-Content -Path $summary -Encoding UTF8
"Cases: d = 45, 50; rng = 1, 2, 3; solvers = SeDuMi, MOSEK" | Add-Content -Path $summary -Encoding UTF8
"Each case runs in a separate MATLAB process." | Add-Content -Path $summary -Encoding UTF8

$ds = @(45, 50)
$seeds = @(1, 2, 3)
$solvers = @("sedumi", "mosek")

foreach ($d in $ds) {
    foreach ($seed in $seeds) {
        foreach ($solver in $solvers) {
            $caseName = "qs_d${d}_rng${seed}_${solver}"
            $logFile = Join-Path $outDir "$caseName.txt"
            $batch = "cd('$workDir'); run_qs_rank_case($d, $seed, '$solver')"

            "[$(Get-Date -Format o)] START $caseName" | Add-Content -Path $summary -Encoding UTF8
            if (Test-Path $logFile) {
                Remove-Item -LiteralPath $logFile -Force
            }

            $matlabArgs = "-wait -logfile `"$logFile`" -batch `"$batch`""
            $proc = Start-Process -FilePath $matlab `
                -ArgumentList $matlabArgs `
                -Wait -PassThru -WindowStyle Hidden

            "[$(Get-Date -Format o)] END $caseName exit=$($proc.ExitCode) log=$logFile" | Add-Content -Path $summary -Encoding UTF8
        }
    }
}

"QS rank sweep launcher finished at $(Get-Date -Format o)" | Add-Content -Path $summary -Encoding UTF8
