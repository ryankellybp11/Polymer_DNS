#!MC 1410

$!VarSet |numDigits| = "6"
$!VarSet |FileBase| = "outputs/flowfield/time-"
$!VarSet |increment| = "10"
$!VarSet |numberoffiles| = "10"
$!VarSet |gridfile| = "/home/rek879/DNS/DNS_2_28_24/outputs/flowfield/grid.szplt"
#$!VarSet |numDigits| = ""
#$!VarSet |FileBase| = ""
#$!VarSet |increment| = ""
#$!VarSet |numberoffiles| = ""

#$!PROMPTFORTEXTSTRING |numDigits|
#    INSTRUCTIONS = "Number of digits in filename"
#$!PROMPTFORTEXTSTRING |FileBase|
#    INSTRUCTIONS = "Base filename"
#$!PROMPTFORTEXTSTRING |increment|
#    INSTRUCTIONS = "Increment of file numbers"
#$!PROMPTFORTEXTSTRING |numberoffiles|
#    INSTRUCTIONS = "Number of total files"


# Import grid and first file
$!ReadDataSet  '"STANDARDSYNTAX" "1.0" "FILELIST_DATAFILES" "2" "|gridfile|" "/home/rek879/DNS/DNS_2_28_24/outputs/flowfield/time-000001.szplt"'
  DataSetReader = 'Tecplot Subzone Data Loader'
  ReadDataOption = New


$!LOOP |numberoffiles|
$!VARSET |add| = (|LOOP|*|increment|)
$!VARSET |n| = "|add|"

# format file numbers to correct digits

$!If |numDigits| == 0
   $!VarSet |finalN| = "|n|"
$!ElseIf |numDigits| == 1
   $!VarSet |finalN| = "|n|" # no use in formatting here
$!ElseIf |numDigits| == 2
  $!VarSet |finalN| = "|n%02d|"
$!ElseIf |numDigits| == 3
  $!VarSet |finalN| = "|n%03d|"
$!ElseIf |numDigits| == 4
  $!VarSet |finalN| = "|n%04d|"
$!ElseIf |numDigits| == 5
  $!VarSet |finalN| = "|n%05d|"
$!ElseIf |numDigits| == 6
  $!VarSet |finalN| = "|n%06d|"
$!Endif

$!PAUSE "|macrofilepath|/|FileBase||finalN|.szplt"

#Checks to see if the file exists
$!EXTENDEDCOMMAND
COMMANDPROCESSORID = "extendmcr"
Command = 'QUERY.FILEEXISTS "|macrofilepath|/|FileBase||finalN|.szplt" "exists"'
$!IF "|exists|" == "YES"
    $!READDATASET '"STANDARDSYNTAX" "1.0" "FILELIST_DATAFILES" "2" "|gridfile|" "|macrofilepath|/|FileBase||finalN|.szplt"'
      DataSetReader = 'Tecplot Subzone Data Loader'
      ReadDataOption = Append
      ResetStyle = No
      AssignStrandIDs = No
      InitialPlotType = Automatic
      InitialPlotFirstZoneOnly = No
      AddZonesToExistingStrands = No
      VarLoadMode = ByName
$!ENDIF

$!ENDLOOP
#$!GlobalRGB RedChannelVar = 4
#$!GlobalRGB GreenChannelVar = 4
#$!GlobalRGB BlueChannelVar = 8
$!SetContourVar 
  Var = 4
  ContourGroup = 1
  LevelInitMode = ResetToNice
#$!SetContourVar 
#  Var = 8
#  ContourGroup = 2
#  LevelInitMode = ResetToNice
#$!SetContourVar 
#  Var = 4
#  ContourGroup = 3
#  LevelInitMode = ResetToNice
#$!SetContourVar 
#  Var = 4
#  ContourGroup = 4
#  LevelInitMode = ResetToNice
#$!SetContourVar 
#  Var = 4
#  ContourGroup = 5
#  LevelInitMode = ResetToNice
#$!SetContourVar 
#  Var = 4
#  ContourGroup = 6
#  LevelInitMode = ResetToNice
#$!SetContourVar 
#  Var = 4
#  ContourGroup = 7
#  LevelInitMode = ResetToNice
$!IsoSurfaceLayers Show = Yes
$!ThreeDView 
  ViewerPosition
    {
    X = 70.3183880278474
    Y = 20.60526973935089
    Z = -22.53995309875651
    }
  ViewWidth = 6.34016
$!IsoSurfaceAttributes 1  Isovalue1 = 25
$!PrintSetup Palette = Color
$!ExportSetup ExportFormat = MPEG4
$!ExportSetup ImageWidth = 653
$!ExportSetup AnimationSpeed = 15
$!ExportSetup ExportFName = './untitled1.mp4'
$!AnimateTime 
  StartTime = 1
  EndTime = 100 
  Skip = 1
  CreateMovieFile = Yes
