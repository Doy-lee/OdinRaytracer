@echo OFF
setlocal
set script_dir=%~dp0
set build_dir=..\Build

mkdir %build_dir% 2>nul
pushd %build_dir%

set compiler=%script_dir%..\..\odin\odin.exe
call %compiler% build %script_dir%\main.odin -debug -show-timings

popd
endlocal
