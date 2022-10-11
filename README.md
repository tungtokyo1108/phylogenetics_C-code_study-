# phylogenetics_C-code_study-

## Install BeagleLib 
### Step 1:
- Follow instruction to download and compile C++ BeagleLib in mac [GitHub BeagleLib](https://github.com/beagle-dev/beagle-lib/wiki/MacInstallInstructions)

### Step 2: Link with program in Xcode
- Choose Xcode > Preferences… from the main menu, then choose the “Locations” tab, then click the “Custom Paths” subtab.
- Click the + button and add a new custom path item named BL_ROOT with display name BL_ROOT. On my computer, I entered the path /usr/local/include/libhmsbeagle-1/
- Click on Build Settings, All (not Basic or Customized), Levels (not Combined), and then type “Header Search Paths” into the Search box to filter out all build settings except the one we are interested in at the moment.
- In the “Header Search Paths” row, double-click on the column with the blue project icon :blueproject: (this is the column for the project settings, the other column labeled strom contains target-specific settings).
- Click the + button at the bottom left corner and type $(BL_ROOT)
- cd ~/lib/static
- cp /usr/local/lib/libhmsbeagle.1.dylib . 
- cp /usr/local/lib/libhmsbeagle.dylib . 
- Click on the Build Phases tab. Click the + sign in Link Binary With Libraries and (after clicking the Add Other… button) navigate to the directory ~/lib/static and select the libhmsbeagle.1.dylib and libhmsbeagle.dylib 
