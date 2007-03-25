@echo off
rem Batch file for timing various pairings
echo A
benchmark < a.param | find "average"
echo D159
benchmark < d159.param | find "average"
echo D201
benchmark < d201.param | find "average"
echo D224
benchmark < d224.param | find "average"
echo E
benchmark < e.param | find "average"
echo F
benchmark < f.param | find "average"
echo G
benchmark < g149.param | find "average"
echo A1
benchmark < a1.param | find "average"
