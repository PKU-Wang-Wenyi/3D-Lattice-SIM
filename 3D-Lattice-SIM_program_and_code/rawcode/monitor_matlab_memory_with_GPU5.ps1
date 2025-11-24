# ===============================
# monitor_matlab_memory_with_GPU_total_v3_fixed_csv_display.ps1
# 功能: 运行 Open_3DSIM.m 并每0.5秒记录 MATLAB CPU 内存 + 私有内存 + GPU 显存
#       实时输出到屏幕并保存为 CSV
# ===============================

# === 日志文件 ===
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logFile = "matlab_memory_gpu_log_$timestamp.csv"
if (Test-Path $logFile) { Remove-Item $logFile }

# === 写入表头 ===
"Timestamp,WorkingSet(GB),PrivateMemory(GB),GPU_Mem_Used(GB)" | Out-File -FilePath $logFile -Encoding utf8

# === 启动 MATLAB ===
$matlab = Start-Process -FilePath "matlab.exe" -ArgumentList "-batch `"run('Main_3ISIM_ZT_AF.m');exit;`"" -PassThru
Write-Host "Start monitoring MATLAB (PID=$($matlab.Id))"
Write-Host "Logging to $logFile`n"

# === 初始化峰值记录 ===
$peakWorkingSetGB = 0
$peakPrivateMemGB = 0
$peakGpuGB = 0

while ($true) {
    $matlabProcs = Get-Process -Name "MATLAB" -ErrorAction SilentlyContinue
    if (-not $matlabProcs) { break }

    try {
        $time = Get-Date -Format "yyyy-MM-dd HH:mm:ss.fff"
        $workingSetGB = [math]::Round(($matlabProcs | Measure-Object -Property WorkingSet64 -Sum).Sum / 1GB, 3)
        $privateMemGB = [math]::Round(($matlabProcs | Measure-Object -Property PrivateMemorySize64 -Sum).Sum / 1GB, 3)

        # GPU 显存
        $gpuUsedGB = 0
        if (Get-Command "nvidia-smi.exe" -ErrorAction SilentlyContinue) {
            $gpuInfo = & nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits
            $gpuUsedGB = [math]::Round((($gpuInfo | ForEach-Object {[int]$_}) | Measure-Object -Sum).Sum / 1024, 3)
        }

        # 格式化输出行
        $line = "$time,$workingSetGB,$privateMemGB,$gpuUsedGB"
        $line | Tee-Object -FilePath $logFile -Append

        # 记录峰值
        if ($workingSetGB -gt $peakWorkingSetGB) { $peakWorkingSetGB = $workingSetGB }
        if ($privateMemGB -gt $peakPrivateMemGB) { $peakPrivateMemGB = $privateMemGB }
        if ($gpuUsedGB -gt $peakGpuGB) { $peakGpuGB = $gpuUsedGB }

    } catch {
        break
    }

    Start-Sleep -Milliseconds 500
}

# === 写入峰值信息 ===
"`nPeak Summary,,," | Tee-Object -FilePath $logFile -Append
"Peak WorkingSet (GB),$peakWorkingSetGB,,," | Tee-Object -FilePath $logFile -Append
"Peak PrivateMemory (GB),$peakPrivateMemGB,,," | Tee-Object -FilePath $logFile -Append
"Peak GPU Memory Used (GB),$peakGpuGB,,," | Tee-Object -FilePath $logFile -Append

Write-Host "`nMonitoring finished. Log saved to $logFile"
