from datetime import datetime
from launchtle import launch_tle


launch_date_ref = datetime(2018, 11, 11, 3, 0)
tle_ref = ["1 70002U 18599A   18315.16151858  .00000000  00000-0  00000-0 0    07",
           "2 70002  85.1205  90.1568 0012705 292.5520 107.9249 15.20792276    04"]

# Testing the Implementation against 'launchtle' from sattools
# launchtle -c irvine.txt -i 70002 -t 2018-11-11T03:00:00 -T 2018-11-11T03:00:00 -I 70002 -d 18599A
fixtures = [(datetime(2018, 11, 11, 3, 0),
            ["1 70002U 18599A   18315.16151858  .00000000  00000-0  00000-0 0    07",
             "2 70002  85.1205  90.1568 0012705 292.5520 107.9249 15.20792276    04"]),
            (datetime(2018, 11, 11, 3, 5),
            ["1 70002U 18599A   18315.16499080  .00000000  00000-0  00000-0 0    09",
             "2 70002  85.1205  91.4102 0012705 292.5520 107.9249 15.20792276    02"]),
            (datetime(2018, 11, 11, 4, 0),
            ["1 70002U 18599A   18315.20318525  .00000000  00000-0  00000-0 0    08",
             "2 70002  85.1205 105.1979 0012705 292.5520 107.9249 15.20792276    07"]),
            (datetime(2018, 11, 12, 3, 0),
            ["1 70002U 18599A   18316.16151858  .00000000  00000-0  00000-0 0    08",
             "2 70002  85.1205  91.1424 0012705 292.5520 107.9249 15.20792276    06"])]

def test_launchtle():
    for new_launch_date, tle_correct in fixtures:
        tle = launch_tle(['DUMMYSAT', *tle_ref], launch_date_ref, new_launch_date)
        assert(tle_correct[1] == tle[2])
        assert(tle_correct[0] == tle[1])
