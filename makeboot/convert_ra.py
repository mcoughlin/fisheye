
import numpy as np

def ra2decimal(hour,minute,second):
    angle_radians = (hour + (minute/60.0) + (second/3600.0))*(2*np.pi)/24.0
    return angle_radians*360.0/(2*np.pi)

def dec2decimal(hour,minute,second):
    if hour < 0:
        return (hour - (minute/60.0) - (second/3600.0))
    else:   
        return (hour + (minute/60.0) + (second/3600.0))

hour_ra = 18.0
minute_ra = 36.0
second_ra = 56.5

hour_dec = 38.0
minute_dec = 47.0
second_dec = 6.2

print ra2decimal(hour_ra,minute_ra,second_ra), dec2decimal(hour_dec,minute_dec,second_dec)

