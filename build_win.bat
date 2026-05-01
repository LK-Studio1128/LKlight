@chcp 65001 >nul
@echo off
REM Build LKlight for Windows (native)
REM ============================================================
REM 编译 LKlight 单二进制版本（v1.0+，含 cmd_generate rayon 并行优化）
REM 输出：target\release\LKlight.exe
REM ============================================================
REM 推荐 GNU (MinGW-w64) 工具链（无需 MSVC）：
REM   方案 A - Scoop:   scoop install mingw
REM   方案 B - MSYS2:   pacman -S mingw-w64-x86_64-gcc
REM   方案 C - WinLibs: https://winlibs.com/  (解压后将 bin 加入 PATH)
REM   方案 D - Conda:   conda install -c conda-forge m2w64-toolchain
REM
REM 安装 MinGW-w64 后执行：
REM   rustup target add x86_64-pc-windows-gnu
REM   然后运行本脚本
REM
REM 如已有 MSVC (Visual Studio Build Tools 2019+)，本脚本会自动回退到 MSVC。
REM ============================================================

setlocal enabledelayedexpansion

cd /d "%~dp0"

echo === LKlight Windows build ===
rustc --version
cargo --version

set BUILD_OK=0
set TARGET_GNU=x86_64-pc-windows-gnu
set BINARY_GNU=%~dp0target\%TARGET_GNU%\release\LKlight.exe
set BINARY_MSVC=%~dp0target\release\LKlight.exe

REM 优先尝试 GNU 工具链（无需 MSVC）
rustup target list --installed 2>nul | findstr /C:"%TARGET_GNU%" >nul
if not errorlevel 1 (
    echo.
    echo [GNU] 使用 x86_64-pc-windows-gnu 工具链编译...
    cargo build --release --target %TARGET_GNU%
    if not errorlevel 1 (
        set BINARY=%BINARY_GNU%
        set BUILD_OK=1
        echo [GNU] 编译成功
    ) else (
        echo [GNU] 编译失败，尝试 MSVC 工具链...
    )
) else (
    echo [GNU] 未安装 GNU 目标，尝试添加...
    rustup target add %TARGET_GNU% 2>nul
    rustup target list --installed 2>nul | findstr /C:"%TARGET_GNU%" >nul
    if not errorlevel 1 (
        cargo build --release --target %TARGET_GNU%
        if not errorlevel 1 (
            set BINARY=%BINARY_GNU%
            set BUILD_OK=1
            echo [GNU] 编译成功
        )
    )
)

REM 回退到 MSVC
if "%BUILD_OK%"=="0" (
    echo.
    echo [MSVC] 尝试 MSVC 默认工具链...
    cargo build --release
    if errorlevel 1 (
        echo ERROR: 所有工具链编译均失败
        echo 请安装 MinGW-w64（推荐）或 Visual Studio Build Tools
        exit /b 1
    )
    set BINARY=%BINARY_MSVC%
    set BUILD_OK=1
    echo [MSVC] 编译成功
)

echo.
echo === Build successful ===
echo Binary: %BINARY%

REM 打包到 dist\LKlight-win\
set DIST=%~dp0dist\LKlight-win
if exist "%DIST%" rmdir /s /q "%DIST%"
mkdir "%DIST%"
copy "%BINARY%" "%DIST%\LKlight.exe"

echo.
echo === Package ready ===
echo Output: %DIST%\
echo   LKlight.exe   (单二进制，含全部子命令；约 17 MB)
echo.
echo 子命令:
echo   setup      ^<rec.pdb^> ^<lig.pdb^> ^<method^> [-s N] [-g N] [--anm]
echo   run        ^<setup.json^> ^<initial_positions.dat^> ^<steps^> ^<method^>
echo   generate   ^<rec.pdb^> ^<lig.pdb^> ^<gso.out^> ^<N^>     (rayon 并行)
echo   cluster    ^<gso.out^>
echo   rank       ^<num_swarms^> ^<steps^>
echo   top        ^<ranking.list^> ^<N^>
echo   filter     ^<ranking.list^> ^<restraints.list^>
echo   gso_to_csv ^<ranking.list^> ^<out.csv^>
echo   score      ^<rec.pdb^> ^<lig.pdb^> ^<method^> [--tx X --ty Y --tz Z]
echo   trajectory ^<rec.pdb^> ^<lig.pdb^> ^<swarm_id^> ^<glowworm_id^> ^<steps^>
echo   pipeline   ^<rec.pdb^> ^<lig.pdb^> ^<method^> [--threads N]
echo.
echo Methods: dfire fastdfire dfire2 dna mj3h pydock cpydock sd vdw pisa sipper tobi ddna
echo.
echo 部署：复制 LKlight.exe 到 LKDock 引擎目录，例如：
echo   ${LKDock_v3.0_Win}\engine\LKlight\LKlight.exe
echo.
echo 快速测试:
echo   "%DIST%\LKlight.exe" --help
echo.

endlocal
