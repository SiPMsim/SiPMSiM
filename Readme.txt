############# References ############# 

This is a standalone version of a SiPM class initially developed to be incorporated into GATE. 
GATE and Geant4 are not required to run this code. 

Implementation: B. Mehadji et al., Monte Carlo simulation of SiPMs with GATE, J. Instrum. 17 (2022) P09025

Experimental validation: B. Mehadji et al., Monte Carlo simulation of a scintillation crystal read by a SiPM with GATE, Nucl. Instrum. Meth. A 1048 (2023) 167905 


############# Compilation ############# 


Required packages:  libxml2-dev libgsl-dev zlib1g-dev
Compile using cmake :
create a build folder of custom name (for instance SiPM_build)
open terminal in the build folder
execute: 

cmake "SiPM code folder directory"  
make 
make install 


############# To run ############# 

Go to build folder and type command

./SiPM "SiPM.xml file location" "name of SiPM in the xml file"  " file location of photons list to be processed in NPY format (should include x(in mm),y(in mm),z (in mm),time (in ns) and optionally energy (in ev)) 


Exemple: 

Copy exemple folder to build folder
Go to build folder
Create a folder named result
open terminal in folder and type command

./SiPM ./exemple  mppc-13360-3050  ./exemple/XYZTime_ev.npy ./result

Run the jupyter notbook in exemple folder to see the result


############# SiPM.xml file instantiation ############# 

###### Optional paramaters #######  

PDE, SEED; SPTR

###### Custom parameters (one should use one of them at ounce. If several are present in SiPM.xml, the first one will be taken into account) #####

#### Signal shape ###

Bi-exponential function exemple: 

<property name="tauRise" value="1" unit="nanosecond"/>
<property name="tauFall" value="100" unit="nanosecond"/>


Spice model from F. Corsi et al., Nucl. Instrum. Methods Phys. Res. A (2007): 

<property name="circuit" Rs="value"  Rq="value"  Cd="value" Cq="value" Cg="value"/>

Sampled signal exemple: 

<propertyvector name="PULSE" unit="nanosecond" >
	<ve time="0.0" value="0.011675983255021347" ></ve>
	<ve time="0.05000000066756627" value="0.014392492543427967" ></ve>
	<ve time="0.10000000133514675" value="0.01854140934527878" ></ve>
	<ve time="0.15000000200271302" value="0.022908229910304254" ></ve>
	...
</propertyvector>


### Optical Crosstalk geometrical distribution ###

Distance based (from firing microcell and sum should be equal to one) exemple:

<propertyvector name="CROSSTALK_DISPERTION" >
	<ve value="9.32349329e-01"></ve>
	<ve value="6.16959618e-02"></ve>
	<ve value="5.38537476e-03"></ve>
	<ve value="5.13045742e-04"></ve>
	<ve value="5.06332074e-05"></ve>
	<ve value="5.09200439e-06"></ve>
	<ve value="5.12085052e-07"></ve>
	<ve value="5.14986007e-08"></ve>
</propertyvector>


Map based (sum equal to one, square of non pair dimensions, centre of square represents the firing microcell and probability should be set to zero) example:

<propertyvector name="CROSSTALK_MAP"  dim="3" >
	<ve X="1" Y="1" value="0.125" ></ve>
	<ve X="1" Y="2" value="0.125" ></ve>
	<ve X="1" Y="3" value="0.125" ></ve>
	<ve X="2" Y="1" value="0.125" ></ve>
	<ve X="2" Y="2" value="0" ></ve>
	<ve X="2" Y="3" value="0.125" ></ve>
	<ve X="3" Y="1" value="0.125" ></ve>
	<ve X="3" Y="3" value="0.125" ></ve>
	<ve X="3" Y="3" value="0.125" ></ve>                                     
</propertyvector>

###### For more details about implementation of noise, refer to ######

B. Mehadji et al., Monte Carlo simulation of SiPMs with GATE, J. Instrum. 17 (2022) P09025

############# Optical photons #############

If compiled as is, NPY file can be read with columns in the fallowing order:

x(in mm),y(in mm),z (in mm),time (in ns) and optionally energy (in ev)


Important ! : 

Optical photons should be within the SiPM, code will exit if is not the case with out of bound error
If one mendatory parameter is not present in SiPM.xml the code will also not run. 
One can set noise values to zero if not intended to be simulated but they still need to be described in SiPM.xml.

ANY CONTRIBUTION IS WELCOMED !
