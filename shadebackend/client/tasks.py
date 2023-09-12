from background_task import background

@background(schedule=60*60*60)  # Schedule the task to run every hour
def my_scheduled_task():
    # Your task code goes here
    print("Executing scheduled task now...")
    # This function will be executed at the scheduled time
