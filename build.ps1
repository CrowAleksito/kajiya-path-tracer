$ErrorActionPreference = "Stop"

# Buscar g++ en PATH o en la instalacion de WinGet LLVM-MinGW
$gpp = Get-Command g++ -ErrorAction SilentlyContinue | Select-Object -ExpandProperty Source
if (-not $gpp) {
    $wingetLinks = "$env:LOCALAPPDATA\Microsoft\WinGet\Links\g++.exe"
    if (Test-Path $wingetLinks) {
        $gpp = $wingetLinks
    } else {
        $found = Get-ChildItem "$env:LOCALAPPDATA\Microsoft\WinGet\Packages" -Recurse -Filter "x86_64-w64-mingw32-g++.exe" -ErrorAction SilentlyContinue | Select-Object -First 1 -ExpandProperty FullName
        $gpp = $found
    }
}

if (-not $gpp) {
    Write-Error "No se encontro g++."
    exit 1
}

Write-Host "Compilando con: $gpp" -ForegroundColor Cyan
New-Item -ItemType Directory -Force -Path build | Out-Null

& $gpp -std=c++17 -O3 -Wall -Wextra -static -Iinclude src/main.cpp -lgdi32 -luser32 -pthread -o build/kajiya_pathtracer.exe
& $gpp -std=c++17 -O3 -Wall -Wextra -static -Iinclude tests/test_math_and_integrals.cpp -pthread -o build/test_math_and_integrals.exe

Write-Host "Compilacion exitosa: build/kajiya_pathtracer.exe y build/test_math_and_integrals.exe" -ForegroundColor Green
