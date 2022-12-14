all: clean compile copy jar

clean: 
	rm -f NGSEPcore_4.3.0.jar
	rm -rf bin
	
compile:
	mkdir bin 
	javac -cp lib/jsci-core.jar:lib/htsjdk-2.22.jar -d bin src/ngsep/*.java src/ngsep/*/*.java src/ngsep/*/*/*.java

copy: 
	cp -f src/ngsep/transcriptome/ProteinTranslatorDefaultBundle.properties bin/ngsep/transcriptome/
	cp -f src/ngsep/main/CommandsDescriptor.xml bin/ngsep/main/
	cp -f src/ngsep/genome/*.js bin/ngsep/genome/

jar: 
	mkdir dist
	jar -xf lib/jsci-core.jar JSci
	mv JSci dist/
	jar -xf lib/htsjdk-2.22.jar htsjdk
	mv htsjdk dist/
	cp -r bin/* dist/
	jar -cfe NGSEPcore_4.3.0.jar ngsep.NGSEPcore -C dist . 
	rm -rf dist
