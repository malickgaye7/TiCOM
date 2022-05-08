## Installation
1. In the terminal, navigate to the directory into which you would like to download TiCOM.
2. Clone the TiCOM repository using the command
    ```
    git clone https://github.com/malickgaye7/TiCOM.git
    ```
3. Run TiCOM in its installation directory, ```../TiCOM/```, using the ```./TiCOM``` command. If you don't know what to do, please run ```./TiCOM help``` or consult the documentation.

    It is recommended (but in no means necessary) to add the TiCOM directory to ```PATH``` such that TiCOM can be easily run without having to reference the TiCOM installation directory (if you don't know how to do this, see this article: https://linuxize.com/post/how-to-add-directory-to-path-in-linux/).
    
## Basic Usage
Once installed, TiCOM can be called using the ```./TiCOM``` command, assuming that the user is in the TiCOM installation directory. If the user is in some other directory (and does not have the TiCOM directory added to their path), then the user can access TiCOM via

    
    [TiCOM installation directory]/TiCOM
    
... otherwise, if the TiCOM directory is in ```PATH``` then a user can call TiCOM by simply using the ```TiCOM``` command.


\
TiCOM comes with two modules: Citcom & TiRADE. TiCOM starts by using a specified initial Citcom input file; an example is provided in the Citcom directory, ```./Citcom_axi_for_TiCOM/input-example-Citcom-JKedit2```. This file contains information about the body being simulated.
   
   * TiRADE & Citcom can also be run independently - example input files for each module is contained in the appropriate module directory.

\
The ```./settings``` file more so specifies how TiCOM and its modules should be run. This file specifies:
   * TiRADE and Citcom's locations (such that these modules can be switched out if needed)
   * Initial Citcom input file
   * The first step number, last step number, and step increment (amount by which the step increases between runs)
   * Data output parameters
   * Various TiRADE parameters

Only variable *assignments* should be edited in this file. Any variable that ends in 'flag' should take on a value of either ```0``` or ```1``` (which corresponds to off/on).

\
The ```./functions``` file contains critical utility functions for TiCOM. It is not advised to make removals from this file as that may affect TiCOM's functionality. This file may be sourced in the terminal independently of running TiCOM.


\
Please see the documentation for more about how to run TiCOM.
