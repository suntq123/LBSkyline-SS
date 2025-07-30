# LBSkyline-SS
This repo is the code for the paper "Beyond Access Pattern Privacy: Efficient Volume-Hiding Location-Based Skyline Queries".

## Requirements
* JDK-17

## Get Started
For ease of testing, we have combined the code for the data user and the two servers into one project. The test requires three networked computers C1,  C2 and a client User.
Before running the program, import the project into C1, C2 and User. All parameter configurations can be completed in the main function of /src/main/java/cn/ac/iscas/TestSSQV1.java. 
      
## Query Processing
* Computer C1
    * C1 selects the role "c1" by modifying the code in the main function of /src/main/java/cn/ac/iscas/TestSSQV1.java as follows:  
        ```
        //select a role  
        args = c1.split(" ");  
        //args = c2.split(" ");  
        //args = user.split(" ");
        ```  
    * C1 runs the TestSSQV1.java.
* Computer C2
    * C2 selects the role "c2" by modifying the code in the main function of /src/main/java/cn/ac/iscas/TestSSQV1.java as follows:  
        ```
        //select a role  
        //args = c1.split(" ");  
        args = c2.split(" ");  
        //args = user.split(" ");  
        ```
    * C2 runs the TestSSQV1.java.
* Computer User
    * User selects the role "user" by modifying the code in the main function of /src/main/java/cn/ac/iscas/TestSSQV1.java as follows: 
        ``` 
        //select a role  
        //args = c1.split(" ");  
        //args = c2.split(" ");  
        args = user.split(" ");  
        ```
    * User runs the TestSSQV1.java.
