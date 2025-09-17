#!/bin/bash
#change the following parameters
runningTimeLimitHours=3 #hours
sleepTime=30            #seconds
#  ========================================================
username=$(whoami)                                 # Set the username for which to check the jobs
runningTime=0                                      # Initialize the running time to 0
condorQFile='condor.out'                           # File to save the output of condor_q
runningTimeLimit=$((runningTimeLimitHours * 3600)) #seconds
flag=1                                             # Flag to check if it is the first time

# Function to check if all jobs for the user are finished
check_jobs() {
    condor_q >$condorQFile # Save the output of condor_q to a file

    # Extract the line containing the user's job information
    userJobsLine=$(grep "Total for $username" $condorQFile)
    # Extract the number of jobs from the user's line
    leftJobs=$(echo "$userJobsLine" | awk '{print $4}')     # 4th column (word) in the line
    idleJobs=$(echo "$userJobsLine" | awk '{print $10}')    # 10th column (word) in the line
    runningJobs=$(echo "$userJobsLine" | awk '{print $12}') # 12th column (word) in the line
    heldJobs=$(echo "$userJobsLine" | awk '{print $14}')    # 14th column (word) in the line

    #example userJobLine:
    #Total for username: 10 jobs; 1 completed, 2 removed, 3 idle, 0 running, 4 held, 0 suspended
    if [ "$flag" = "1" ]; then
        totalJobs=leftJobs
        flag=0
    fi

    echo "$userJobsLine"
    # Check if all jobs are finished
    if [ "$leftJobs" = "0" ]; then
        echo "Jobs for user $username are finished! -------------------> DONE"
        return 1
    # Check if the running time exceeds the limit
    elif [ $runningTime -gt $runningTimeLimit ]; then
        echo "Jobs for user $username are running for too long! Removing all remaining jobs..."
        condor_rm $username
        echo "Further resubmission needed"
        return 1
    # Check if there are held jobs
    elif [ "$heldJobs" != "0" ] && [ "$idleJobs" = "0" ] && [ "$runningJobs" = "0" ]; then
        echo "There are held jobs for user $username. Removing all remaining jobs..."
        condor_rm $username
        echo "Further resubmission needed"
        return 1

    else
        # If there are still running jobs, wait for 10 seconds
        echo "Some jobs are still running for user $username. Please wait."
        echo "Sleep $sleepTime seconds..."
        if [ "$runningJobs" != "0" ]; then
            runningTime=$((runningTime + sleepTime)) # Increment runningTime by 10 seconds
        fi
        # output percentage of time passed
        echo "Percentage of time passed: $((runningTime * 100 / runningTimeLimit))%"
        echo "Percentage of jobs done: $((100 - leftJobs * 100 / totalJobs))%"
        runningTimeHours=$((runningTime / 3600))
        echo "Total running time: $runningTimeHours hours"
        echo " "
        sleep $sleepTime
    fi
}
# Main loop to continuously check job status
while :; do
    check_jobs || break # If all jobs are finished, exit the loop
done

rm $condorQFile
