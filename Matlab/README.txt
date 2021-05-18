# AusRecon

Code Repository of AusRecon, a Prior Austenite reconstruction code for written
in MATLAB  by Alexander Brust, Steve Niezgoda, and Eric Payton. This repository
is maintained by Austin Gerlt, who is not an author, just an enthusiastic user.

NOTES ON GETTING DATA FOR THIS CODE:
While any .ang file can (in theory) be used,  a large set of useful DISTRO A 
datasets can be found at either of the locations listed below:
AF Research Lab link: 
https://drive.google.com/drive/folders/1cw21aKhuaXX_yEq_0JA0GJTm7L0Oi91o?usp=sharing
ELSZ Storage File location:
sftp://sshfs.rdte.afrl.dren.mil/   location: /project/RXCM/DISTRO_A_Datasets_AusRecon

NOTES ON THE CODE ITSELF:

1)  This code was primarily written by Alexander Brust as part of his PhD thesis, 
    with significant contributions from Eric Payton at AFRL and Seve Niezgoda at
    Ohio State University. It uses a graph cutting algorithm to assign martensitic 
    Block, Packet, and Variant IDs to every pixel scanned on a given EBSD map. For 
    Further questions on the specifics of how this algorithm has been applied, 
    contact Alex Brust (alexander.brust.ctr@us.af.mil).

2)  the reconstruction code is still in mid-beta, and as such has some bugs.
    Most notably, there are occasional dropped pixels (approximately 5 per scan)
	Users should also note that Martensitic structures associated with twinned 
	grains have different IDs in the code's output than non-twinned grains.

3)  If users utilize the 347 DISTRO A EBSD scans mentioned above, it should be
	known that these scans were collected by Vikas Sinha (vikas.sinha.1.ctr@us.af.mil)
	and Jared Shank (Jared.Shank.1.ctr@us.af.mil). **Please credit them as well 
	as everyone else mentioned above appropriately when using this data.** 
