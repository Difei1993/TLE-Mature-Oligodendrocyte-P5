###############################################################################
# 32B-1
# Balanced cNMF pilot using H5AD inputs
#
# Steps:
#   prepare -> factorize -> combine -> k_selection_plot
#
# Consensus is NOT run.
###############################################################################

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

###############################################################################
# 0. Paths
###############################################################################

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

$WorkRoot = Join-Path `
    $Project `
    "04_tables\32B_CNMf_WORK_20260728_215801"

$H5adRoot = Join-Path `
    $WorkRoot `
    "_h5ad_inputs"

$LogRoot = Join-Path `
    $WorkRoot `
    "_logs_h5ad"

$ReviewRoot = Join-Path `
    $WorkRoot `
    "_K_SELECTION_REVIEW_H5AD"

$OverallLog = Join-Path `
    $LogRoot `
    "32B1_overall_run.log"

###############################################################################
# 1. Environment audit
###############################################################################

foreach ($RequiredPath in @(
    $Project,
    $EnvPath,
    $PythonExe,
    $CnmfExe,
    $H5adRoot
)) {
    if (-not (Test-Path -LiteralPath $RequiredPath)) {
        throw "Required path not found: $RequiredPath"
    }
}

New-Item `
    -ItemType Directory `
    -Force `
    -Path $WorkRoot | Out-Null

New-Item `
    -ItemType Directory `
    -Force `
    -Path $LogRoot | Out-Null

"32B-1 H5AD cNMF pilot" |
    Set-Content `
        -Path $OverallLog `
        -Encoding UTF8

"Started: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" |
    Add-Content `
        -Path $OverallLog `
        -Encoding UTF8

Write-Host ""
Write-Host "============================================================"
Write-Host "32B-1 cNMF pilot using H5AD inputs"
Write-Host "============================================================"

& $PythonExe --version

if ($LASTEXITCODE -ne 0) {
    throw "Python environment validation failed."
}

# Validate cNMF without truncating its help-output pipe.
$CnmfVersion = & $PythonExe -c "import importlib.metadata as m; print(m.version('cnmf'))"

if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($CnmfVersion)) {
    throw "cNMF package validation failed."
}

Write-Host "cNMF version: $CnmfVersion"

###############################################################################
# 2. Lineage definitions
###############################################################################

$LineagePlan = @(
    @{
        Name = "Astrocytes"
        InputFile = Join-Path $H5adRoot "Astrocytes_primary.h5ad"
        KValues = @(4..12)
    },

    @{
        Name = "Microglia_macrophage"
        InputFile = Join-Path $H5adRoot "Microglia_macrophage_primary.h5ad"
        KValues = @(4..12)
    },

    @{
        Name = "OPC"
        InputFile = Join-Path $H5adRoot "OPC_primary.h5ad"
        KValues = @(3..10)
    },

    @{
        Name = "Oligodendrocytes"
        InputFile = Join-Path $H5adRoot "Oligodendrocytes_primary.h5ad"
        KValues = @(4..14)
    }
)

###############################################################################
# 3. H5AD input validation
###############################################################################

function Test-H5adInput {

    param(
        [Parameter(Mandatory = $true)]
        [string]$Lineage,

        [Parameter(Mandatory = $true)]
        [string]$InputFile
    )

    if (-not (
        Test-Path `
            -LiteralPath $InputFile `
            -PathType Leaf
    )) {
        throw "$Lineage H5AD file not found: $InputFile"
    }

    $InputItem = Get-Item `
        -LiteralPath $InputFile

    if ($InputItem.Length -le 0) {
        throw "$Lineage H5AD file is empty: $InputFile"
    }

    Write-Host ""
    Write-Host "$Lineage input:"
    Write-Host "  File: $InputFile"
    Write-Host "  Size: $([math]::Round($InputItem.Length / 1MB, 2)) MB"
}

Write-Host ""
Write-Host "============================================================"
Write-Host "Preflight H5AD validation"
Write-Host "============================================================"

foreach ($LineageSpec in $LineagePlan) {
    Test-H5adInput `
        -Lineage $LineageSpec.Name `
        -InputFile $LineageSpec.InputFile
}

Write-Host ""
Write-Host "All four H5AD inputs passed validation."

###############################################################################
# 4. Robust native-command runner
#
# Start-Process is used so that Python warnings written to stderr are not
# incorrectly converted into terminating PowerShell NativeCommandError events.
###############################################################################

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

    foreach ($TempFile in @(
        $StdoutTemp,
        $StderrTemp
    )) {
        if (Test-Path $TempFile) {
            Remove-Item `
                -LiteralPath $TempFile `
                -Force
        }
    }

    Write-Host ""
    Write-Host "============================================================"
    Write-Host $Label
    Write-Host "Started: $($StartTime.ToString('yyyy-MM-dd HH:mm:ss'))"
    Write-Host "============================================================"

    "Command: $CnmfExe $($Arguments -join ' ')" |
        Set-Content `
            -Path $LogFile `
            -Encoding UTF8

    $Process = Start-Process `
        -FilePath $CnmfExe `
        -ArgumentList $Arguments `
        -Wait `
        -PassThru `
        -NoNewWindow `
        -RedirectStandardOutput $StdoutTemp `
        -RedirectStandardError $StderrTemp

    if (Test-Path $StdoutTemp) {

       $StdoutLines = @(
    Get-Content `
        -LiteralPath $StdoutTemp `
        -ErrorAction SilentlyContinue
)

if (@($StdoutLines).Count -gt 0) {

    $StdoutLines |
        Add-Content `
            -Path $LogFile `
            -Encoding UTF8

    $StdoutLines |
        ForEach-Object {
            Write-Host $_
        }
}
    }

    if (Test-Path $StderrTemp) {

$StderrLines = @(
    Get-Content `
        -LiteralPath $StderrTemp `
        -ErrorAction SilentlyContinue
)

if (@($StderrLines).Count -gt 0) {

    "" |
        Add-Content `
            -Path $LogFile `
            -Encoding UTF8

    "===== STDERR =====" |
        Add-Content `
            -Path $LogFile `
            -Encoding UTF8

    $StderrLines |
        Add-Content `
            -Path $LogFile `
            -Encoding UTF8

    $StderrLines |
        ForEach-Object {
            Write-Host $_
        }
}
    }

    $EndTime = Get-Date
    $Elapsed = New-TimeSpan `
        -Start $StartTime `
        -End $EndTime

    if ($Process.ExitCode -ne 0) {

        "FAILED: $Label; exit code=$($Process.ExitCode)" |
            Add-Content `
                -Path $OverallLog `
                -Encoding UTF8

        throw @"
$Label failed.
Exit code: $($Process.ExitCode)
Log file: $LogFile
"@
    }

    "COMPLETED: $Label; elapsed=$($Elapsed.ToString())" |
        Add-Content `
            -Path $OverallLog `
            -Encoding UTF8

    foreach ($TempFile in @(
        $StdoutTemp,
        $StderrTemp
    )) {
        if (Test-Path $TempFile) {
            Remove-Item `
                -LiteralPath $TempFile `
                -Force
        }
    }

    Write-Host ""
    Write-Host "Completed: $Label"
    Write-Host "Elapsed: $($Elapsed.ToString())"
}

###############################################################################
# 5. Run one lineage
###############################################################################

function Run-CnmfLineage {

    param(
        [Parameter(Mandatory = $true)]
        [string]$Lineage,

        [Parameter(Mandatory = $true)]
        [string]$InputFile,

        [Parameter(Mandatory = $true)]
        [int[]]$KValues
    )

    $OutputDirectory = Join-Path `
        $WorkRoot `
        $Lineage

    $RunName = "${Lineage}_primary"

    New-Item `
        -ItemType Directory `
        -Force `
        -Path $OutputDirectory | Out-Null

    $PrepareMarker = Join-Path `
        $OutputDirectory `
        "32B1_PREPARE_COMPLETED.txt"

    $FactorizeMarker = Join-Path `
        $OutputDirectory `
        "32B1_FACTORIZE_COMPLETED.txt"

    $CombineMarker = Join-Path `
        $OutputDirectory `
        "32B1_COMBINE_COMPLETED.txt"

    $KSelectionMarker = Join-Path `
        $OutputDirectory `
        "32B1_K_SELECTION_COMPLETED.txt"

    $PrepareLog = Join-Path `
        $LogRoot `
        "${Lineage}_01_prepare.log"

    $FactorizeLog = Join-Path `
        $LogRoot `
        "${Lineage}_02_factorize.log"

    $CombineLog = Join-Path `
        $LogRoot `
        "${Lineage}_03_combine.log"

    $KSelectionLog = Join-Path `
        $LogRoot `
        "${Lineage}_04_k_selection.log"

    $KArguments = @(
        $KValues |
        ForEach-Object {
            [string]$_
        }
    )

    Write-Host ""
    Write-Host "############################################################"
    Write-Host "Lineage: $Lineage"
    Write-Host "K values: $($KValues -join ', ')"
    Write-Host "Input: $InputFile"
    Write-Host "############################################################"

    ###########################################################################
    # prepare
    ###########################################################################

    if (-not (Test-Path $PrepareMarker)) {

        $PrepareArguments = @(
            "prepare",
            "--output-dir", $OutputDirectory,
            "--name", $RunName,
            "-c", $InputFile,
            "-k"
        ) + $KArguments + @(
            "--n-iter", "30",
            "--total-workers", "1",
            "--seed", "20260728",
            "--numgenes", "2000",
            "--max-nmf-iter", "2000"
        )

        Invoke-CnmfStep `
            -Label "$Lineage : prepare" `
            -LogFile $PrepareLog `
            -Arguments $PrepareArguments

        "prepare completed" |
            Set-Content `
                -Path $PrepareMarker `
                -Encoding ASCII

    } else {
        Write-Host "$Lineage : prepare already completed; skipping."
    }

    ###########################################################################
    # factorize
    ###########################################################################

    if (-not (Test-Path $FactorizeMarker)) {

        $FactorizeArguments = @(
            "factorize",
            "--output-dir", $OutputDirectory,
            "--name", $RunName,
            "--worker-index", "0",
            "--total-workers", "1",
            "--skip-completed-runs"
        )

        Invoke-CnmfStep `
            -Label "$Lineage : factorize" `
            -LogFile $FactorizeLog `
            -Arguments $FactorizeArguments

        "factorize completed" |
            Set-Content `
                -Path $FactorizeMarker `
                -Encoding ASCII

    } else {
        Write-Host "$Lineage : factorize already completed; skipping."
    }

    ###########################################################################
    # combine
    ###########################################################################

    $CombineArguments = @(
        "combine",
        "--output-dir", $OutputDirectory,
        "--name", $RunName
    )

    Invoke-CnmfStep `
        -Label "$Lineage : combine" `
        -LogFile $CombineLog `
        -Arguments $CombineArguments

    "combine completed" |
        Set-Content `
            -Path $CombineMarker `
            -Encoding ASCII

    ###########################################################################
    # k-selection plot
    ###########################################################################

    $KSelectionArguments = @(
        "k_selection_plot",
        "--output-dir", $OutputDirectory,
        "--name", $RunName
    )

    Invoke-CnmfStep `
        -Label "$Lineage : k_selection_plot" `
        -LogFile $KSelectionLog `
        -Arguments $KSelectionArguments

    "k selection completed" |
        Set-Content `
            -Path $KSelectionMarker `
            -Encoding ASCII

$KFiles = @(
    Get-ChildItem `
        -Path $OutputDirectory `
        -Recurse `
        -File `
        -ErrorAction SilentlyContinue |
    Where-Object {
        $_.Name -match "k_selection"
    }
)

    if ($KFiles.Count -eq 0) {
        throw "$Lineage completed but no K-selection output was found."
    }
}

###############################################################################
# 6. Run four lineages
###############################################################################

$OverallStart = Get-Date

foreach ($LineageSpec in $LineagePlan) {

    Run-CnmfLineage `
        -Lineage $LineageSpec.Name `
        -InputFile $LineageSpec.InputFile `
        -KValues $LineageSpec.KValues
}

$OverallEnd = Get-Date

$OverallElapsed = New-TimeSpan `
    -Start $OverallStart `
    -End $OverallEnd

###############################################################################
# 7. Build review directory
###############################################################################

if (Test-Path $ReviewRoot) {
    Remove-Item `
        -LiteralPath $ReviewRoot `
        -Recurse `
        -Force
}

New-Item `
    -ItemType Directory `
    -Force `
    -Path $ReviewRoot | Out-Null

foreach ($LineageSpec in $LineagePlan) {

    $Lineage = $LineageSpec.Name

    $LineageWork = Join-Path `
        $WorkRoot `
        $Lineage

    $LineageReview = Join-Path `
        $ReviewRoot `
        $Lineage

    New-Item `
        -ItemType Directory `
        -Force `
        -Path $LineageReview | Out-Null

    $KSelectionFiles = Get-ChildItem `
        -Path $LineageWork `
        -Recurse `
        -File `
        -ErrorAction SilentlyContinue |
    Where-Object {
        $_.Name -match "k_selection"
    }

    foreach ($File in $KSelectionFiles) {
        Copy-Item `
            -LiteralPath $File.FullName `
            -Destination (
                Join-Path `
                    $LineageReview `
                    $File.Name
            ) `
            -Force
    }

    $LineageLogs = Get-ChildItem `
        -Path $LogRoot `
        -File `
        -Filter "${Lineage}_*.log" `
        -ErrorAction SilentlyContinue

    foreach ($LogFile in $LineageLogs) {
        Copy-Item `
            -LiteralPath $LogFile.FullName `
            -Destination (
                Join-Path `
                    $LineageReview `
                    $LogFile.Name
            ) `
            -Force
    }
}

Copy-Item `
    -LiteralPath $OverallLog `
    -Destination (
        Join-Path `
            $ReviewRoot `
            "32B1_overall_run.log"
    ) `
    -Force

###############################################################################
# 8. Review summary
###############################################################################

$ReviewSummary = Join-Path `
    $ReviewRoot `
    "32B1_K_SELECTION_REVIEW_SUMMARY.txt"

@"
32B-1 H5AD cNMF K-selection pilot completed.

Started:
$($OverallStart.ToString("yyyy-MM-dd HH:mm:ss"))

Completed:
$($OverallEnd.ToString("yyyy-MM-dd HH:mm:ss"))

Elapsed:
$($OverallElapsed.ToString())

Lineages:
- Astrocytes: K 4-12
- Microglia/macrophage: K 4-12
- OPC: K 3-10
- Oligodendrocytes: K 4-14

NMF replicates:
30 per K.

High-variance genes:
2000.

Discovery restrictions:
- No senescence signature used.
- No disease label used.
- Cells were balanced across all eight donors.

Consensus:
NOT RUN.

Next step:
Review stability and reconstruction error before selecting K.
"@ |
    Set-Content `
        -Path $ReviewSummary `
        -Encoding UTF8

###############################################################################
# 9. ZIP and SHA256
###############################################################################

$ReviewZip = Join-Path `
    $WorkRoot `
    "32B1_K_SELECTION_REVIEW_H5AD.zip"

if (Test-Path $ReviewZip) {
    Remove-Item `
        -LiteralPath $ReviewZip `
        -Force
}

Compress-Archive `
    -Path (
        Join-Path `
            $ReviewRoot `
            "*"
    ) `
    -DestinationPath $ReviewZip `
    -CompressionLevel Optimal `
    -Force

$Hash = Get-FileHash `
    -LiteralPath $ReviewZip `
    -Algorithm SHA256

$HashValue = $Hash.Hash.ToLower()

$HashLine = "$HashValue  $([IO.Path]::GetFileName($ReviewZip))"

$HashFile = "${ReviewZip}.sha256.txt"

Set-Content `
    -Path $HashFile `
    -Value $HashLine `
    -Encoding ASCII

###############################################################################
# 10. Final marker
###############################################################################

$CompletionFile = Join-Path `
    $WorkRoot `
    "32B1_H5AD_COMPLETED_SUCCESSFULLY.txt"

@"
32B-1 H5AD pilot completed successfully.

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
        -Path $CompletionFile `
        -Encoding UTF8

Write-Host ""
Write-Host "============================================================"
Write-Host "32B-1 H5AD cNMF pilot completed successfully"
Write-Host "============================================================"
Write-Host "Elapsed:"
Write-Host $OverallElapsed
Write-Host ""
Write-Host "K-selection review ZIP:"
Write-Host $ReviewZip
Write-Host ""
Write-Host "SHA256:"
Write-Host $HashLine
Write-Host ""
Write-Host "Consensus was NOT run."
Write-Host "============================================================"