Repository for handling scripts and scratches for GEOFRAME/OMS3 project
=======================================================================

**Next steps for the project are**:
- [x] Create a script to run compilation, project directory structure generation and simfiles
- [x] Test the script thorougfully
- [] Think about the implementation of some FORTRAN code, starting from the links provided:
	1. https://alm.engr.colostate.edu/cb/wiki/16993
	2. https://alm.engr.colostate.edu/cb/wiki/17097
	3. https://alm.engr.colostate.edu/cb/wiki/17105#section-OMS3_2FFortran9x+Integration+Example_
	4. https://teamwork.niwa.co.nz/display/IFM/Implementing+a+FORTRAN+model+in+OMS3 > see sparrow example
- [x] Think about the needed JAVA classes for the Richards integration code from *Casulli et al*.
- [x] Test OMS3 and R integration 
- [] Port the Gumbel code in an OMS3 module through R

**Should I take a look at some FORTRAN components that we already want to integrate inside the GEOFRAME?**
- [HydroGen] (https://github.com/geoframecomponents/Hydro_gen), a computer code for generating two-dimensional space random functions with an assigned covariance structure. The original code is written in Ansi Fortran 77.
- ...

Links
=====
* [Markdown syntax] (https://guides.github.com/features/mastering-markdown/)
* [GEOFrame source code](https://github.com/geogramecomponents)
* [GEOFrame projects](https://github.com/GEOframeOMSProjects)
* [Last Renjin doclink seen](http://docs.renjin.org/en/latest/library/moving-data-between-java-and-r-code.html)

Trivia
======
* This has been a pain... **REMEMBER** to include the base class path inside the `-classpath` during JAVA program execution! *e.g.*:
`java -classpath .:../lib/renjin-script-engine-0.8.2165-jar-with-dependencies.jar TryRenjin`
`.` obviously stands for "this dir", so the commant above assumes you're executing the code inside the class directory. 
* **THE BIGGEST PAINSOURCE OF ALL**: for some reason, `java` does not interpret correctly `~` as `$HOME`. So, better use `$HOME` directly for scripting. 
* Example of packaged program execution (also, nice use of dot notation for upwards directory navigation): `java -classpath ./:../../../JAVA_BASE_CLASSES/renjin-script-engine-0.8.2165-jar-with-dependencies.jar example1.TryRenjin` 

Remember: *five pushes a day keeps the doctor away*


Integration test between Java-OMS3-R
===================================

This has been only a playground for testing the integration of OMS3 with R
statistical language through an R wrapper. 

**Chosen backends**:
* Definitely going to use *[Renjin](http://www.renjin.org/)*. Seems promising:

> "The biggest advantage of Renjin is that the R interpreter itself is a Java module which can be seamlessly integrated into any Java application. This dispenses with the need to load dynamic libraries or to provide some form of communication between separate processes.These types of interfaces are often the source of much agony because they place very specific demands on the environment in which they run."

* [TUTORIAL](http://docs.renjin.org/en/latest/introduction.html)
* [INSTALLATION PROCEDURE]()
	* `sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys B71EB5FA4FE5A5F6`
	* `sudo su`
	* `echo "deb https://nexus.bedatadriven.com/content/repositories/renjin-release/ ./" > /etc/apt/sources.list.d/renjin.list`
	* `apt-get update`
	* `apt-get install renjin`
* [PACKAGES](http://docs.renjin.org/en/latest/interactive/index.html): 
> From within Renjin’s REPL, there is no `install.packages()` function: the first time you try to load a package with `library()`, Renjin will check the repository for a package with the matching name and download it to a local repository located in `~/.m2/repository`.

* From [PROJECT SETUP](http://docs.renjin.org/en/latest/library/project-setup.html), it seems that for maximum flexibility I'll have to use a single JAR file if I'm not willing to use any project organizer (like Maven, or Eclipse). ANT should suffice as a build tool, with the automatically generated project directory //TODO: verify this
* [EXAMPLE CODE](http://docs.renjin.org/en/latest/library/evaluating.html)
		  
> From within Renjin’s REPL, there is no install.packages() function: the first time you try to load a package with library(), Renjin will check the repository for a package with the matching name and download it to a local repository located in ~/.m2/repository.


**Tested backends**:
* *[JRI](https://rforge.net/JRI/)*, which works through [REngine](https://github.com/s-u/REngine).
	* No compilation problem reported, but class does not seem to instantiate correctly
	during class launch. Also, lacks of thoroug documentation.
	
	Base class to be loaded:
	```java
	import org.rosuda.JRI.Rengine;
	```
**Backends in testing**:
* -

**PROGRESS - PROJECT COMPLETED**
- [x] Tested R and JAVA integration succesfully;
- [x] Test R and OMS3 integration by writing one OMS component which takes as input a strimg from the _simfile_
