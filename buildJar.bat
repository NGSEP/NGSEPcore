mkdir dist
jar -xf lib/jsci-core.jar JSci
move JSci dist/
jar -xf lib/htsjdk-2.22.jar htsjdk
move htsjdk dist/
xcopy /E bin dist
jar -cfe NGSEPcore_5.1.1.jar ngsep.NGSEPcore -C dist . 
rmdir /s /q dist
:END
cmd /k