mkdir dist
jar -xf lib/jsci-core.jar JSci
move JSci dist/
jar -xf lib/htsjdk-1.129.jar htsjdk
move htsjdk dist/
xcopy /E bin dist
jar -cfe NGSEPcore_3.3.3.jar ngsep.NGSEPcore -C dist . 
rmdir /s /q dist
:END
cmd /k