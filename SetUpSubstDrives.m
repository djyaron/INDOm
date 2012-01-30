%% Unload substituted drives

system('subst s: /d');
system('subst o: /d');
system('subst p: /d');

% If "Invalid parameter" error, then the drives didn't exist as subst drives.

%% Load substituted drives

rootdir = 'C:\Users\Christian\Documents\Research\Yaron\dyes2\data\5-merPPV\';
SFolder = 'Exp';
OFolder = 'INDOLib';
PFolder = 'GSLib';

system(['subst s: "',rootdir,SFolder,'"']);
system(['subst o: "',rootdir,OFolder,'"']);
system(['subst p: "',rootdir,PFolder,'"']);