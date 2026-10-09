###############################################################################
# 32B-2A
# Independent balanced-sampling confirmation of candidate K values
#
# No consensus.
# No disease labels or senescence signatures used.
###############################################################################

$ErrorActionPreference = "Stop"

$Project = $env:P5_WORK_ROOT

$EnvPath = Join-Path `
    $Project `
    "04_tables\python_envs\cnmf310"

$PythonExe = Join-Path `
    $EnvPath `
    "Scripts\python.exe"

$CnmfExe = Join-Path `
    $EnvPath `
    "Scripts\cnmf.exe"

$PrimaryWorkRoot = Join-Path `
    $Project `
    "04_tables\32B_CNMf_WORK_20260728_215801"

$WorkRoot = Join-Path `
    $Project `
    "04_tables\32B2A_CNMf_SENSITIVITY_20260728_215801"

$H5adRoot = Join-Path `
    $WorkRoot `
    "_h5ad_inputs"

$LogRoot = Join-Path `
    $WorkRoot `
    "_logs"

$ReviewRoot = Join-Path `
    $WorkRoot `
    "_K_CONFIRMATION_REVIEW"

$OverallLog = Join-Path `
    $LogRoot `
    "32B2A_overall_run.log"

foreach ($RequiredPath in @(
    $PythonExe,
    $CnmfExe,
    $H5adRoot,
    $PrimaryWorkRoot
)) {
    if (-not (Test-Path -LiteralPath $RequiredPath)) {
        throw "Required path not found: $RequiredPath"
    }
}

New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
New-Item -ItemType Directory -Force -Path $LogRoot | Out-Null

$Version = & $PythonExe -c `
    "import importlib.metadata as m; print(m.version('cnmf'))"

if (
    $LASTEXITCODE -ne 0 -or
    [string]::IsNullOrWhiteSpace($Version)
) {
    throw "cNMF package validation failed."
}

Write-Host "Python:"
& $PythonExe --version
Write-Host "cNMF version: $Version"

"32B-2A independent K confirmation" |
    Set-Content -Path $OverallLog -Encoding UTF8

"Started: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" |
    Add-Content -Path $OverallLog -Encoding UTF8


$LineagePlan = @(
    @{
        Name = "Astrocytes"
        Input = Join-Path $H5adRoot "Astrocytes_sensitivity.h5ad"
        K = @(7, 8)
    },
    @{
        Name = "Microglia_macrophage"
        Input = Join-Path $H5adRoot "Microglia_macrophage_sensitivity.h5ad"
        K = @(6, 7)
    },
    @{
        Name = "OPC"
        Input = Join-Path $H5adRoot "OPC_sensitivity.h5ad"
        K = @(6, 7)
    },
    @{
        Name = "Oligodendrocytes"
        Input = Join-Path $H5adRoot "Oligodendrocytes_sensitivity.h5ad"
        K = @(6, 7)
    }
)


function Invoke-CnmfStep {

    param(
        [Parameter(Mandatory = $true)]
        [string]$Label,

        [Parameter(Mandatory = $true)]
        [string]$LogFile,

        [Parameter(Mandatory = $true)]
        [string[]]$Arguments
    )

    $StartTime = Get-Date

    $StdoutTemp = "${LogFile}.stdout.tmp"
    $StderrTemp = "${LogFile}.stderr.tmp"

    foreach ($TempFile in @($StdoutTemp, $StderrTemp)) {
        if (Test-Path $TempFile) {
            Remove-Item -LiteralPath $TempFile -Force
        }
    }

    Write-Host ""
    Write-Host "============================================================"
    Write-Host $Label
    Write-Host "Started: $($StartTime.ToString('yyyy-MM-dd HH:mm:ss'))"
    Write-Host "============================================================"

    "Command: $CnmfExe $($Arguments -join ' ')" |
        Set-Content -Path $LogFile -Encoding UTF8

    $Process = Start-Process `
        -FilePath $CnmfExe `
        -ArgumentList $Arguments `
        -Wait `
        -PassThru `
        -NoNewWindow `
        -RedirectStandardOutput $StdoutTemp `
        -RedirectStandardError $StderrTemp

    $StdoutLines = @(
        Get-Content `
            -LiteralPath $StdoutTemp `
            -ErrorAction SilentlyContinue
    )

    $StderrLines = @(
        Get-Content `
            -LiteralPath $StderrTemp `
            -ErrorAction SilentlyContinue
    )

    if (@($StdoutLines).Count -gt 0) {
        $StdoutLines |
            Add-Content -Path $LogFile -Encoding UTF8

        $StdoutLines |
            ForEach-Object { Write-Host $_ }
    }

    if (@($StderrLines).Count -gt 0) {
        "" | Add-Content -Path $LogFile -Encoding UTF8
        "===== STDERR =====" |
            Add-Content -Path $LogFile -Encoding UTF8

        $StderrLines |
            Add-Content -Path $LogFile -Encoding UTF8

        $StderrLines |
            ForEach-Object { Write-Host $_ }
    }

    if ($Process.ExitCode -ne 0) {
        throw "$Label failed; exit code=$($Process.ExitCode); log=$LogFile"
    }

    $EndTime = Get-Date
    $Elapsed = New-TimeSpan -Start $StartTime -End $EndTime

    "COMPLETED: $Label; elapsed=$($Elapsed.ToString())" |
        Add-Content -Path $OverallLog -Encoding UTF8

    foreach ($TempFile in @($StdoutTemp, $StderrTemp)) {
        if (Test-Path $TempFile) {
            Remove-Item -LiteralPath $TempFile -Force
        }
    }

    Write-Host "Completed: $Label"
    Write-Host "Elapsed: $($Elapsed.ToString())"
}


function Run-Lineage {

    param(
        [string]$Lineage,
        [string]$InputFile,
        [int[]]$KValues
    )

    if (-not (
        Test-Path -LiteralPath $InputFile -PathType Leaf
    )) {
        throw "$Lineage input not found: $InputFile"
    }

    $OutputDirectory = Join-Path $WorkRoot $Lineage
    $RunName = "${Lineage}_sensitivity"

    New-Item `
        -ItemType Directory `
        -Force `
        -Path $OutputDirectory | Out-Null

    $PrepareMarker = Join-Path `
        $OutputDirectory `
        "32B2A_PREPARE_COMPLETED.txt"

    $FactorizeMarker = Join-Path `
        $OutputDirectory `
        "32B2A_FACTORIZE_COMPLETED.txt"

    $KArguments = @(
        $KValues |
        ForEach-Object { [string]$_ }
    )

    if (-not (Test-Path $PrepareMarker)) {
        $Arguments = @(
            "prepare",
            "--output-dir", $OutputDirectory,
            "--name", $RunName,
            "-c", $InputFile,
            "-k"
        ) + $KArguments + @(
            "--n-iter", "30",
            "--total-workers", "1",
            "--seed", "20260729",
            "--numgenes", "2000",
            "--max-nmf-iter", "2000"
        )

        Invoke-CnmfStep `
            -Label "$Lineage : prepare" `
            -LogFile (
                Join-Path $LogRoot "${Lineage}_01_prepare.log"
            ) `
            -Arguments $Arguments

        "completed" |
            Set-Content -Path $PrepareMarker -Encoding ASCII
    }

    if (-not (Test-Path $FactorizeMarker)) {
        $Arguments = @(
            "factorize",
            "--output-dir", $OutputDirectory,
            "--name", $RunName,
            "--worker-index", "0",
            "--total-workers", "1",
            "--skip-completed-runs"
        )

        Invoke-CnmfStep `
            -Label "$Lineage : factorize" `
            -LogFile (
                Join-Path $LogRoot "${Lineage}_02_factorize.log"
            ) `
            -Arguments $Arguments

        "completed" |
            Set-Content -Path $FactorizeMarker -Encoding ASCII
    }

    Invoke-CnmfStep `
        -Label "$Lineage : combine" `
        -LogFile (
            Join-Path $LogRoot "${Lineage}_03_combine.log"
        ) `
        -Arguments @(
            "combine",
            "--output-dir", $OutputDirectory,
            "--name", $RunName
        )

    Invoke-CnmfStep `
        -Label "$Lineage : k_selection_plot" `
        -LogFile (
            Join-Path $LogRoot "${Lineage}_04_k_selection.log"
        ) `
        -Arguments @(
            "k_selection_plot",
            "--output-dir", $OutputDirectory,
            "--name", $RunName
        )

    $KFiles = @(
        Get-ChildItem `
            -Path $OutputDirectory `
            -Recurse `
            -File |
        Where-Object {
            $_.Name -match "k_selection"
        }
    )

    if (@($KFiles).Count -eq 0) {
        throw "$Lineage has no K-selection output."
    }
}


$OverallStart = Get-Date

foreach ($Spec in $LineagePlan) {
    Run-Lineage `
        -Lineage $Spec.Name `
        -InputFile $Spec.Input `
        -KValues $Spec.K
}

$OverallEnd = Get-Date
$OverallElapsed = New-TimeSpan `
    -Start $OverallStart `
    -End $OverallEnd


if (Test-Path $ReviewRoot) {
    Remove-Item `
        -LiteralPath $ReviewRoot `
        -Recurse `
        -Force
}

New-Item -ItemType Directory -Force -Path $ReviewRoot | Out-Null

foreach ($Spec in $LineagePlan) {
    $Lineage = $Spec.Name

    $LineageReview = Join-Path $ReviewRoot $Lineage
    $PrimaryReview = Join-Path $LineageReview "primary"
    $SensitivityReview = Join-Path $LineageReview "sensitivity"

    New-Item -ItemType Directory -Force -Path $PrimaryReview | Out-Null
    New-Item -ItemType Directory -Force -Path $SensitivityReview | Out-Null

    $PrimaryDirectory = Join-Path $PrimaryWorkRoot $Lineage
    $SensitivityDirectory = Join-Path $WorkRoot $Lineage

    $PrimaryFiles = @(
        Get-ChildItem `
            -Path $PrimaryDirectory `
            -Recurse `
            -File `
            -ErrorAction SilentlyContinue |
        Where-Object {
            $_.Name -match "k_selection"
        }
    )

    foreach ($File in $PrimaryFiles) {
        Copy-Item `
            -LiteralPath $File.FullName `
            -Destination (
                Join-Path $PrimaryReview $File.Name
            ) `
            -Force
    }

    $SensitivityFiles = @(
        Get-ChildItem `
            -Path $SensitivityDirectory `
            -Recurse `
            -File `
            -ErrorAction SilentlyContinue |
        Where-Object {
            $_.Name -match "k_selection"
        }
    )

    foreach ($File in $SensitivityFiles) {
        Copy-Item `
            -LiteralPath $File.FullName `
            -Destination (
                Join-Path $SensitivityReview $File.Name
            ) `
            -Force
    }

    $Logs = @(
        Get-ChildItem `
            -Path $LogRoot `
            -File `
            -Filter "${Lineage}_*.log" `
            -ErrorAction SilentlyContinue
    )

    foreach ($LogFile in $Logs) {
        Copy-Item `
            -LiteralPath $LogFile.FullName `
            -Destination (
                Join-Path $SensitivityReview $LogFile.Name
            ) `
            -Force
    }
}


$SummaryPath = Join-Path `
    $ReviewRoot `
    "32B2A_K_CONFIRMATION_SUMMARY.txt"

@"
32B-2A independent balanced-sampling K confirmation completed.

Elapsed:
$($OverallElapsed.ToString())

Candidate K values:
- Astrocytes: 7, 8
- Microglia/macrophage: 6, 7
- OPC: 6, 7
- Oligodendrocytes: 6, 7

Replicates:
30 per K.

Sensitivity sampling seed:
20260729.

Discovery restrictions:
- No disease label used.
- No senescence signature used.
- Eight donors equally represented.

Consensus:
NOT RUN.

Next:
Select one final K per lineage only after comparing primary and
independent-sensitivity stability and reconstruction error.
"@ |
    Set-Content -Path $SummaryPath -Encoding UTF8


$ReviewZip = Join-Path `
    $WorkRoot `
    "32B2A_K_CONFIRMATION_REVIEW.zip"

if (Test-Path $ReviewZip) {
    Remove-Item -LiteralPath $ReviewZip -Force
}

Compress-Archive `
    -Path (Join-Path $ReviewRoot "*") `
    -DestinationPath $ReviewZip `
    -CompressionLevel Optimal `
    -Force

$Hash = Get-FileHash `
    -LiteralPath $ReviewZip `
    -Algorithm SHA256

$HashLine = "$($Hash.Hash.ToLower())  $([IO.Path]::GetFileName($ReviewZip))"

Set-Content `
    -Path "${ReviewZip}.sha256.txt" `
    -Value $HashLine `
    -Encoding ASCII

@"
32B-2A completed successfully.

Completed:
$($OverallEnd.ToString("yyyy-MM-dd HH:mm:ss"))

Elapsed:
$($OverallElapsed.ToString())

Review ZIP:
$ReviewZip

SHA256:
$HashLine

Consensus was not run.
"@ |
    Set-Content `
        -Path (
            Join-Path $WorkRoot "32B2A_COMPLETED_SUCCESSFULLY.txt"
        ) `
        -Encoding UTF8

Write-Host ""
Write-Host "============================================================"
Write-Host "32B-2A completed successfully"
Write-Host "============================================================"
Write-Host "Elapsed: $OverallElapsed"
Write-Host ""
Write-Host "Review ZIP:"
Write-Host $ReviewZip
Write-Host ""
Write-Host "SHA256:"
Write-Host $HashLine
Write-Host ""
Write-Host "Consensus was NOT run."
Write-Host "============================================================"