mkdir dist
jar -xf lib/jsci-core.jar JSci
move JSci dist/
jar -xf lib/htsjdk-1.129.jar htsjdk
move htsjdk dist/
xcopy /E bin dist
jar -cfe NGSEPcore_4.0.2.jar ngsep.NGSEPcore -C dist . 
rmdir /s /q dist
:END
cmd /k