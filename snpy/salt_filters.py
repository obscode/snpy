'''A database of SALT filters and how to convert them to SNooPy filters.'''

snpy_to_salt = {
   # CSPI
   'B':('@SWOPE2','B','@VEGAHST'),
   'V':('@SWOPE2','V','@VEGAHST'),
   'u':('@SWOPE2','u','@BD17-CSP'),
   'g':('@SWOPE2','g','@BD17-CSP'),
   'r':('@SWOPE2','r','@BD17-CSP'),
   'i':('@SWOPE2','i','@BD17-CSP'),
   # CSPII
   'B2':('@SWOPE3','B','@VEGAHST'),
   'V2':('@SWOPE3','V','@VEGAHST'),
   'u2':('@SWOPE3','u','@BD17-CSP'),
   'g2':('@SWOPE3','g','@BD17-CSP'),
   'r2':('@SWOPE3','r','@BD17-CSP'),
   'i2':('@SWOPE3','i','@BD17-CSP'),
   # Standard
   'Us':('@STANDARD','U','@VEGAHST'),
   'Bs':('@STANDARD','B','@VEGAHST'),
   'Vs':('@STANDARD','V','@VEGAHST'),
   'Rs':('@STANDARD','R','@VEGAHST'),
   'Is':('@STANDARD','I','@VEGAHST'),
   # SDSS
   'u_s':('@SDSS','u','@BD17-CSP'),
   'g_s':('@SDSS','g','@BD17-CSP'),
   'r_s':('@SDSS','r','@BD17-CSP'),
   'i_s':('@SDSS','i','@BD17-CSP'),
   'z_s':('@SDSS','z','@BD17-CSP'),
   # KEPLERCAM
   'Bk':('@KEPLERCAM','B','@VEGAHST'),
   'Vk':('@KEPLERCAM','V','@VEGAHST'),
   'uk':('@KEPLERCAM','u','@BD17-CSP'),
   'rk':('@KEPLERCAM','r','@BD17-CSP'),
   'ik':('@KEPLERCAM','i','@BD17-CSP'),
   # 4-Shooter
   'U4sh':('@4SHOOTER2','U','@VEGAHST'),
   'B4sh':('@4SHOOTER2','B','@VEGAHST'),
   'V4sh':('@4SHOOTER2','V','@VEGAHST'),
   'R4sh':('@4SHOOTER2','R','@VEGAHST'),
   'I4sh':('@4SHOOTER2','U','@VEGAHST'),
   #USNO
   'u_40':('@USNO','u','@BD17-CSP'),
   'g_40':('@USNO','g','@BD17-CSP'),
   'r_40':('@USNO','r','@BD17-CSP'),
   'i_40':('@USNO','i','@BD17-CSP'),
   'z_40':('@USNO','z','@BD17-CSP'),
}   
   
# Stock SALT2 config
snpy_to_salt0 = snpy_to_salt.copy()
for filt in ['u','g','r','i','B','V']:
   snpy_to_salt0[filt] = ('@SWOPE',filt,'@VEGA2')

