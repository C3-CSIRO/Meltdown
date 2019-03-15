IF EXIST C:\Anaconda\python.exe (
C:\Anaconda\python.exe %~dp0\source\MeltdownBatch.py
) ELSE (
C:\Anaconda2\python.exe %~dp0\source\MeltdownBatch.py
)
pause