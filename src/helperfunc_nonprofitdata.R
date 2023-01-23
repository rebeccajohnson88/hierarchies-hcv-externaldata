



clean_select_nccs <- function(data, which_type){
  
  ## filter out outside of scope outside of ors
  ## and select vars of interest
  df_touse = data %>% filter(OUTNCCS == "IN") %>% #apply recommended filter to in-scope
    mutate(which_type = which_type) %>%
    dplyr::select(EIN, CENSUSTRACT, 
                  NTEECC)
  return(df_touse)
  
}


#The following function will download and prepare NCCS Core Files for analysis.
#The function requires the following fields:
#     year: year of data (as of this writing, available years include 1989-2015 for most data sets)
#     type: type of data, given as follows (include quotes)
#           "pc" = 501(c)(3) public charities
#           "pf" = 501(c)(3) private foundations
#           "co" = all other 501(c) organizations

# From here: https://nccs.urban.org/sites/default/files/2018-10/Prep%20NCCS%20Core%20File_1.R
#Create function to download core files
getcorefile <- function(coreyear, coretype) {
  #ensure coretype is lowercase
  coretype <- tolower(coretype)
  
  #make sure user entered pc, pf, or co
  check <- isTRUE(coretype %in% c("pc", "pf", "co"))
  if (check == "FALSE"){
    stop("Please make sure to enter \"pc\" (for public charities), \"pf\" (for private foundations)', or \"co\" (for other 501(c) organizations)")
  }
  
  #create URL for download based on year of file
  if(coretype =="co"){
    URL <- paste("http://nccs-data.urban.org/data/core/", as.character(coreyear), "/coreco.core", as.character(coreyear),as.character(coretype), ".csv", sep ="") 
  }
  if(coretype  %in% c("pc", "pf")) {
    URL <- paste("http://nccs-data.urban.org/data/core/", as.character(coreyear), "/nccs.core", as.character(coreyear),as.character(coretype), ".csv", sep ="")
  }
  
  #get raw data
  rawcore <- GET(as.character(URL))
  if (rawcore$status_code != 200){
    stop("Please check database: it is possible that Core year/type combination you entered does not exist.  See http://nccs-data.urban.org/data.php?ds=core for list of all available datasets, and enter a valid year/type combination")
  }
  stop_for_status(rawcore)
  
  #read data as CSV file format, with proper formatting
  
  #Note 1: the code below reads in most fields from the relevant Core file.  However, users should consider commenting out lines
  #that they do not require.  By default, fields used in NCCS validation processes are not included in output file.
  #Note 2: The code below assumes column names for most recent years of the Core files.  Older editions of the NCCS Core Files
  #may have other column names.  Please consult the relevant data dictionaries (at http://nccs-data.urban.org/data-dictionaries.php)
  #for more information, and change as neecssary.  
  
  
  #parse PC and CO files
  if(coretype  %in% c("pc", "co")) {
    cfile <- content(rawcore, type = "text/csv", 
                     col_types = cols_only(  EIN = col_character(),
                                             ACCPER = col_character(),
                                             ACTIV1 = col_character(),
                                             ACTIV2 = col_character(),
                                             ACTIV3 = col_character(),
                                             ADDRESS = col_character(),
                                             AFCD = col_character(),
                                             ASS_BOY = col_double(),
                                             ASS_EOY = col_double(),
                                             #BLOCK = col_character(),
                                             BOND_BOY = col_double(),
                                             BOND_EOY = col_double(),
                                             CENSUSTRACT = col_character(),
                                             CITY = col_character(),
                                             CLASSCD = col_character(),
                                             COMPENS = col_double(),
                                             #COMPENSP = col_double(),
                                             CONT = col_double(),
                                             CONTACT = col_character(),
                                             #CONTP = col_double(),
                                             DEDUCTCD = col_character(),
                                             DIREXP = col_double(),
                                             EOSTATUS = col_character(),
                                             EPOSTCARD = col_character(),
                                             EXPS = col_double(),
                                             #EXPSP = col_double(),
                                             FILENAME = col_character(),
                                             FIPS = col_character(),
                                             FISYR = col_double(),
                                             #FISYRP = col_double(),
                                             FNDNCD = col_character(),
                                             FRCD = col_character(),
                                             FUNDBAL = col_double(),
                                             FUNDFEES = col_double(),
                                             #FUNDFEESP = col_double(),
                                             FUNDINC = col_double(),
                                             FUNDINCP = col_double(),
                                             GFTGRNTSRCVD170 = col_double(),
                                             GOODS = col_double(),
                                             GRPROF = col_double(),
                                             #GRPROFP = col_double(),
                                             GRREC = col_double(),
                                             GRRECMEM = col_double(),
                                             GRRECOTH = col_double(),
                                             GRRECPUB = col_double(),
                                             GRSINC170 = col_double(),
                                             GRSINCFNDRSNG = col_double(),
                                             GRSINCGAMING = col_double(),
                                             GRSRCPTSADMISSN509 = col_double(),
                                             GRSRCPTSRELATED170 = col_double(),
                                             GRSRNTSPRSNL = col_double(),
                                             GRSRNTSREAL = col_double(),
                                             INITFEES = col_double(),
                                             #INPRIOR = col_integer(),
                                             #INPRIORSRC = col_character(),
                                             INVENTG = col_double(),
                                             INVINC = col_double(),
                                             #INVINCP = col_double(),
                                             LATITUDE = col_double(),
                                             LESSDIRFNDRSNG = col_double(),
                                             LESSDIRGAMING = col_double(),
                                             LEVEL1 = col_character(),
                                             LEVEL2 = col_character(),
                                             LEVEL3 = col_character(),
                                             LEVEL4 = col_character(),
                                             LIAB_BOY = col_double(),
                                             LIAB_EOY = col_double(),
                                             LONGITUDE = col_double(),
                                             MAJGRPB = col_character(),
                                             #MANUALLY_FIXED = col_character(),
                                             MRTG_BOY = col_double(),
                                             MRTG_EOY = col_double(),
                                             MSA_NECH = col_character(),
                                             NAICS = col_character(),
                                             NAME = col_character(),
                                             #NCCSKEY = col_character(),
                                             #NCCSKEY2 = col_character(),
                                             NETA_BOY = col_double(),
                                             NETGNLS = col_double(),
                                             NETINC = col_double(),
                                             NETINCFNDRSNG = col_double(),
                                             NETINCGAMING = col_double(),
                                             #NETINCP = col_double(),
                                             NETRENT = col_double(),
                                             #NETRENTP = col_double(),
                                             NTEE1 = col_character(),
                                             NTEECC = col_character(),
                                             NTEECONF = col_character(),
                                             NTEEFINAL = col_character(),
                                             NTEEFINAL1 = col_character(),
                                             NTEEIRS = col_character(),
                                             NTEESRC = col_character(),
                                             NTMAJ10 = col_character(),
                                             NTMAJ12 = col_character(),
                                             NTMAJ5 = col_character(),
                                             ORGCD = col_character(),
                                             OTHCHGS = col_double(),
                                             OTHINC = col_double(),
                                             OTHINCP = col_double(),
                                             OTHSAL = col_double(),
                                             #OTHSALP = col_double(),
                                             OUTNCCS = col_character(),
                                             OutNCCS = col_character(),
                                             OUTREAS = col_character(),
                                             PAYTAX = col_double(),
                                             PAYTAXP = col_double(),
                                             PMSA = col_character(),
                                             PROGREV = col_double(),
                                             #PROGREVP = col_double(),
                                             Q78A = col_character(),
                                             RANDNUM = col_double(),
                                             REASON = col_character(),
                                             RENTEXP = col_double(),
                                             RENTINC = col_double(),
                                             RETEARN = col_double(),
                                             RETE_BOY = col_double(),
                                             RNTLEXPNSPRSNL = col_double(),
                                             RNTLEXPNSREAL = col_double(),
                                             RNTLINCPRSNL = col_double(),
                                             RNTLINCREAL = col_double(),
                                             ROYALTSINC = col_double(),
                                             RULEDATE = col_character(),
                                             SALEOTHE = col_double(),
                                             SALEOTHG = col_double(),
                                             SALEOTHN = col_double(),
                                             #SALEOTHNP = col_double(),
                                             SALESECN = col_double(),
                                             #SALESECNP = col_double(),
                                             SALESEXP = col_double(),
                                             SECUR = col_double(),
                                             SEC_NAME = col_character(),
                                             SOIYR = col_double(),
                                             SOURCE = col_character(),
                                             SPEVTG = col_double(),
                                             SRVCSVAL170 = col_double(),
                                             SRVCSVAL509 = col_double(),
                                             STATE = col_character(),
                                             STYEAR = col_double(),
                                             SUBCD = col_character(),
                                             SUBSECCD = col_character(),
                                             SUBTOTSUPPINC509 = col_double(),
                                             TAXPER = col_integer(),
                                             #TAXPERP = col_double(),
                                             TOTGFTGRNTRCVD509 = col_double(),
                                             TOTREV = col_double(),
                                             TOTREV2 = col_double(),
                                             #TOTREVP = col_double(),
                                             TOTSUPP509 = col_double(),
                                             TXEXMPTBNDSPROCEEDS = col_double(),
                                             TXREVNUELEVIED170 = col_double(),
                                             TXREVNUELEVIED509 = col_double(),
                                             UNSECUREDNOTESEND = col_double(),
                                             #VALIDATION_STATE = col_double(),
                                             ZIP = col_character(),
                                             ZIP5 = col_character()))
  }
  
  #parse PF file
  if(coretype == "pf"){
    cfile <- content(rawcore, type = "text/csv", 
                     col_types = cols_only(
                       EIN = col_character(),
                       ACCPER = col_character(),
                       ACTIV1 = col_character(),
                       ACTIV2 = col_character(),
                       ACTIV3 = col_character(),
                       ADDRESS = col_character(),
                       AFCD = col_character(),
                       #ASS_CODE = col_character(),
                       BLOCK = col_character(),
                       CENSUSTRACT = col_character(),
                       CITY = col_character(),
                       CLASSCD = col_character(),
                       CONTACT = col_character(),
                       #DEDUCTCD = col_character(),
                       #EOSTATUS = col_integer(),
                       FIPS = col_character(),
                       FISYR = col_integer(),
                       FNDNCD = col_character(),
                       FRCD = col_character(),
                       #INC_CODE = col_character(),
                       LATITUDE = col_double(),
                       LEVEL1 = col_character(),
                       LEVEL2 = col_character(),
                       LEVEL3 = col_character(),
                       LEVEL4 = col_character(),
                       LONGITUDE = col_double(),
                       MAJGRPB = col_character(),
                       #MANUALLY_FIXED = col_character(),
                       MSA_NECH = col_character(),
                       NAICS = col_character(),
                       NAME = col_character(),
                       #NCCSKEY = col_double(),
                       #NCCSKEY2 = col_double(),
                       NTEE1 = col_character(),
                       NTEECC = col_character(),
                       NTEECONF = col_character(),
                       NTEEFINAL = col_character(),
                       NTEEFINAL1 = col_character(),
                       NTEEIRS = col_character(),
                       NTEESRC = col_character(),
                       NTMAJ10 = col_character(),
                       NTMAJ12 = col_character(),
                       NTMAJ5 = col_character(),
                       ORGCD = col_character(),
                       OUTNCCS = col_character(),
                       OutNCCS = col_character(),
                       OUTREAS = col_character(),
                       P10MIR = col_double(),
                       P10TASST = col_double(),
                       P11DISTR = col_double(),
                       P13UNDST = col_double(),
                       P14A4942 = col_double(),
                       P14ASVLA = col_double(),
                       P14ASVLB = col_double(),
                       P14ASVLC = col_double(),
                       P14ASVLD = col_double(),
                       P14B4942 = col_double(),
                       P14C4942 = col_double(),
                       P14D4942 = col_double(),
                       P14ENDWA = col_double(),
                       P14ENDWB = col_double(),
                       P14ENDWC = col_double(),
                       P14ENDWD = col_double(),
                       P14GINVA = col_double(),
                       P14GINVB = col_double(),
                       P14GINVC = col_double(),
                       P14GINVD = col_double(),
                       P14NADJA = col_double(),
                       P14NADJB = col_double(),
                       P14NADJC = col_double(),
                       P14NADJD = col_double(),
                       P14PSUPA = col_double(),
                       P14PSUPB = col_double(),
                       P14PSUPC = col_double(),
                       P14PSUPD = col_double(),
                       P14QDISA = col_double(),
                       P14QDISB = col_double(),
                       P14QDISC = col_double(),
                       P14QDISD = col_double(),
                       P14T4942 = col_double(),
                       P14TASVL = col_double(),
                       P14TENDW = col_double(),
                       P14TGINV = col_double(),
                       P14TNADJ = col_double(),
                       P14TPSUP = col_double(),
                       P14TQDIS = col_double(),
                       P14TSUPA = col_double(),
                       P14TSUPB = col_double(),
                       P14TSUPC = col_double(),
                       P14TSUPD = col_double(),
                       P14TTSUP = col_double(),
                       P1ADMEXP = col_double(),
                       P1CHAEXP = col_double(),
                       P1CONTPD = col_double(),
                       P1DIVID = col_double(),
                       P1EXCREV = col_double(),
                       P1GINVPF = col_double(),
                       P1GOODS = col_double(),
                       P1GRENTS = col_double(),
                       P1INTREV = col_double(),
                       P1INVEXP = col_double(),
                       P1INVRV = col_double(),
                       P1NADINC = col_double(),
                       P1NADJRV = col_double(),
                       P1NETINV = col_double(),
                       P1NGASTS = col_double(),
                       P1OFCOMP = col_double(),
                       P1OTHINC = col_double(),
                       P1TADJEX = col_double(),
                       P1TCONT = col_double(),
                       P1TEXMEX = col_double(),
                       P1TINVEX = col_double(),
                       P1TOTEXP = col_double(),
                       P1TOTREV = col_double(),
                       P2CASHBV = col_double(),
                       P2CRPBND = col_double(),
                       P2CRPSTK = col_double(),
                       P2EYASST = col_double(),
                       P2EYOTIN = col_double(),
                       P2GVTINV = col_double(),
                       P2MORTG = col_double(),
                       P2MTGNTS = col_double(),
                       P2TASFMV = col_double(),
                       P2TINVSC = col_double(),
                       P2TLIABL = col_double(),
                       P2TOTAST = col_double(),
                       P3EYTFND = col_double(),
                       P6DOMORG = col_character(),
                       P6ESTTX = col_double(),
                       P6EXMPF = col_character(),
                       P6TEXCTX = col_double(),
                       P6TX4940 = col_double(),
                       P6TX511 = col_double(),
                       P6TXCRDT = col_double(),
                       P6TXINV = col_double(),
                       P6TXPNLT = col_double(),
                       P6TXRFD = col_double(),
                       P6TXSUBA = col_double(),
                       P6TXWERR = col_double(),
                       P6TXWTH = col_double(),
                       P7BUSINT = col_character(),
                       P7DISBOR = col_character(),
                       P7DISPAY = col_character(),
                       P7DISSAL = col_character(),
                       P7DISSER = col_character(),
                       P7DISTRA = col_character(),
                       P7ELECT = col_character(),
                       P7INDGRT = col_character(),
                       P7JEOINV = col_character(),
                       P7LIQUID = col_character(),
                       P7ORGGRT = col_character(),
                       P7OTHPUR = col_character(),
                       P7PAYGVT = col_character(),
                       P7POFCLM = col_character(),
                       P7POLIT = col_character(),
                       P7PROPAG = col_character(),
                       P7UNDINC = col_character(),
                       PMSA = col_character(),
                       RANDNUM = col_double(),
                       #RECCODE = col_character(),
                       RULEDATE = col_character(),
                       SEC_NAME = col_character(),
                       SOURCE = col_character(),
                       STATE = col_character(),
                       STYEAR = col_character(),
                       SUBSECCD = col_character(),
                       TAXPER = col_character(),
                       #TFLD = col_integer(),
                       #VALIDATION_STATE = col_double(),
                       ZIP = col_character(),
                       ZIP5 = col_character()
                     ))
  }
  
  
  
  #convert all variable names to uppercase
  names(cfile) <- toupper(names(cfile))
  
  #write output to local drive
  #Note: Users should STRONGLY consider saving the data locally to avoid repeated downloads of the data
  write.csv(cfile, file = as.character(paste("core", as.character(coreyear), as.character(coretype), ".csv", sep="")))
  
  #return output to R for immediate exploration
  return(cfile)
}
