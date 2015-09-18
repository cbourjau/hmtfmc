import os

basepath = "/alice/cern.ch/user/p/pwgpp_mc/2015/17_Week/TestMultiplicity/Test4/quicktest/"
remotefiles = os.popen("alien_find {} galice.root".format(basepath))
print remotefiles
max_files = 1000

for idx, remotefile in enumerate(remotefiles):
    if idx > max_files:
        break
    # break at first empty line
    remotepath = "alien://" + remotefile.strip().replace("galice.root", "")
    if not remotefile:
        break
    try:
        os.mkdir(remotepath.split('/')[-2])
    except OSError:
        print "Not re-downloading remotefile. (Folder already exists)"
        continue
    print os.popen("alien_cp -m -s alien://{} ./{}".format(remotepath+"*.root", remotepath.split('/')[-2])).read()
