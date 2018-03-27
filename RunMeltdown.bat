IF EXIST C:\Anaconda2\python.exe (
	C:\Anaconda2\python.exe source\Meltdown.py
) ELSE IF EXIST %HOMEPATH%\AppData\Local\Continuum\Anaconda2\python.exe (
	%HOMEPATH%\AppData\Local\Continuum\Anaconda2\python.exe source\Meltdown.py
) ELSE IF EXIST C:\Anaconda\python.exe (
	C:\Anaconda\python.exe source\Meltdown.py
) ELSE IF EXIST %HOMEPATH%\AppData\Local\Continuum\Anaconda\python.exe (
	%HOMEPATH%\AppData\Local\Continuum\Anaconda\python.exe source\Meltdown.py
) ELSE (
	python source\Meltdown.py
)
pause