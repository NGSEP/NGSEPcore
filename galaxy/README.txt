NGSEP - Next Generation Sequencing Experience Platform
Version 3.3.1
===========================================================================
Galaxy directory:
=================
Here are the scripts for including NGSEP in a Galaxy environment. Download 
all the files to get the full functionality of NGSEP in your Galaxy instance.

Install NGSEP in a Galaxy instance in five simple steps:
 1. Create a directory called "ngsep" inside the "tools" directory of your 
    Galaxy instance. 
 2. Drop there all the xml scripts from this site. 
 3. Drop there too, the latest jar file from NGSEP: NGSEPcore_3.3.1.jar
 4. Include the following lines in the "config/tool_conf.xml" file.
    
          <section id="ngsep" name="NGS: NGSEP">
             <label id="detection" text="Variants discovery"/>
              <tool file="ngsep/FindVariants.xml"/>
              <tool file="ngsep/MergeVariants.xml"/>
              <tool file="ngsep/MergeVCF.xml"/>
             <label id="vcf_utils" text="VCF manipulation"/>
              <tool file="ngsep/Annotate.xml"/>
              <tool file="ngsep/FilterVCF.xml"/>
              <tool file="ngsep/ConvertVCF.xml"/>
              <tool file="ngsep/ImputeVCF.xml"/>
             <label id="statistics" text="Statistics: BAM/VCF"/>
              <tool file="ngsep/QualStats.xml"/>
              <tool file="ngsep/CoverageStats.xml"/>
              <tool file="ngsep/SummaryStats.xml"/>
              <tool file="ngsep/DiversityStats.xml"/>
              <tool file="ngsep/CompareVCF.xml"/>
             <label id="other" text="Other"/>
              <tool file="ngsep/Generate_sequence_name.xml"/>
              <tool file="ngsep/CompareRD.xml"/>
           </section>

  5. Restart Galaxy and you're ready to enjoy all the new functionalities
     of NGSEP in your own instance.
