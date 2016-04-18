@echo off
rem get some local settings from the registry
rem see http://ss64.com/nt/date.html
rem see http://ss64.com/nt/time.html
for /f "tokens=3" %%d in ('reg query ^"HKEY_CURRENT_USER\Control Panel\International^" /v iCountry ^| find ^"REG_SZ^"') do (set country=%%d)
for /f "tokens=3" %%d in ('reg query ^"HKEY_CURRENT_USER\Control Panel\International^" /v sDate ^| find ^"REG_SZ^"') do (set date_sep=%%d)
for /f "tokens=3" %%d in ('reg query ^"HKEY_CURRENT_USER\Control Panel\International^" /v sTime ^| find ^"REG_SZ^"') do (set time_sep=%%d)
rem get the date in the format YYYY.MM.DD
if %country% equ 1 (
    for /f "tokens=1-5 delims=%date_sep% " %%d in ("%date%") do (set build_date=%%g.%%e.%%f)
) else (
    for /f "tokens=1-5 delims=%date_sep% " %%d in ("%date%") do (set build_date=%%e.%%f.%%g)
)
rem get the time in the format HH:MM
for /f "tokens=1-5 delims=%time_sep% " %%d in ("%time%") do (set build_time=%%d:%%e)
for /f "tokens=1-5 delims=%time_sep%" %%d in ("%time%") do (set build_time_with_leading_space=%%d:%%e)
if /i "%build_time_with_leading_space%" neq "%build_time%" (set build_time=0%build_time%)
rem print build timestamp as used as alternative for PROJECT_RELEASE
echo %build_date% (%build_time%)
