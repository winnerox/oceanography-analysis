@echo off
echo Removing MATLAB right-click menu...
reg delete "HKEY_CLASSES_ROOT\.m\shell\RunMATLAB" /f 2>nul
reg delete "HKEY_CLASSES_ROOT\.m\shell\RunMATLABSilent" /f 2>nul
reg delete "HKEY_CLASSES_ROOT\matlab.m.1\shell\RunMATLAB" /f 2>nul
reg delete "HKEY_CLASSES_ROOT\matlab.m.1\shell\RunMATLABSilent" /f 2>nul
reg delete "HKEY_CLASSES_ROOT\SystemFileAssociations\.m\shell\RunMATLAB" /f 2>nul
reg delete "HKEY_CLASSES_ROOT\SystemFileAssociations\.m\shell\RunMATLABSilent" /f 2>nul
echo Done!
echo Restarting explorer...
taskkill /f /im explorer.exe >nul 2>&1
start explorer.exe
echo Fix completed. Please check if right-click menu is normal.
