import schedule as sc       # https://www.geeksforgeeks.org/python-schedule-library/
import time
import solweig_run

## MARK - IGNORE THIS FILE FOR NOW SINCE WE PLAN TO USE AWS FOR JOBS
def jobs_for_24hrs():
    # run function once
    solweig_run.run_solweig()
    sc.every().hour.do(solweig_run.run_solweig())


sc.every().day.at("18:00").do(jobs_for_24hrs)

while True:

    sc.run_pending()
    time.sleep(1)