## Step 1:pip install and run the software(CAUM)
- `pip install -r requirements.txt`
- `python CAUM.py`
## Step 2:Calculate the chemical ages
- Open xlsx file in main window
- Enter necessary parameters (errors for UO2, ThO2, PbO)
- Click Run to to calculate the chemical ages and errors of individual analytical points
- Save calculated errors and ages to new xlsx file
## Step 3:Analyze ages
- Launch submodule from main window
- Open xlsx file saved above
- Adjust settings as needed
- Click Run to show results
note:Certain functionalities (e.g., `ctypes.windll`) are not supported on Linux or other non-Windows platforms. Users running the software on non-Windows systems may need to comment out or modify the relevant code sections.
