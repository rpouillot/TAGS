library(xtable)



shinyUI(
  #Application title
  fluidPage(
    titlePanel(
          HTML('<meta name="Author" content="Pouillot Regis">
                <meta name="KeyWords" content="validation test absence gold standard test de reference">
                <font color="#000080"><tt>TAGS</tt>:
                 evaluation of </font><font color="#008080">T</font><font color="#000080">ests
               in the </font><font color="#008080">A</font><font color="#000080">bsence
               of a </font><font color="#008080">G</font><font color="#000080">old </font><font color="#008080">S</font><font color="#000080">tandard&nbsp;
               <hr>
                 </font>'),
      windowTitle='Tags: validation of Tests in the Absence of a Gold Standard'
    ),

    HTML('The TAGS software can be used to estimate the sensitivity of two or more diagnostic tests in the absence of a gold standard, provided two or more populations with differing prevalences can be cross-classified based on diagnostic test results. 
          The algorithm in the TAGS software follows the frequentist paradigm and utilizes Newton-Raphson and EM algorithms to generate maximum likelihood estimates.
          The TAGS software can accomodate not only data of the "2 independent tests, 2 populations"-type, but also higher order combinations of numbers of tests and numbers of populations. 
         Recognizing that in some instances, true prevalence may be known for some populations, TAGS is capable of utilizing "reference population data" (where one or more populations is of known disease status).  Parameter estimation using TAGS becomes possible once the number of degrees of freedom given by the data is greater than the number of parameters to be estimated.  A goodness-of-fit test and residual correlations, both of which are provided by TAGS output, provide a means of evaluating model adequacy.
         The algorithm includes 2 strong assumptions: (i) diagnostic tests are assumed to be conditionally independent, and (ii) test diagnostic values are considered constant when applied to different populations.
         <br><br>'),
    HTML("Reference: <a href='http://www.sciencedirect.com/science/article/pii/S0167587701002720'>R. Pouillot, G. Gerbier, I.A. Gardner. ''TAGS'', a program for the evaluation of test accuracy in the absence of a gold standard, Preventive Veterinary Medicine, 53 (1-2), pp 67-81.</a>"),
    
    

hr(),
##############################    
# Formulaire
##############################    

selectInput("qtest", HTML("<font color='#000080'><b>1: Enter the number of test(s) to be evaluated</b></font>"), choices = 0:10, selected=2),
selectInput("qpop", HTML("<font color='#000080'><b>2: Enter the number of tested population(s) with an unknown infection status</b></font>"), choices = 0:10, selected=2),
selectInput("qref", HTML("<font color='#000080'><b>3: Enter the number and the category of the reference population(s) tested</b></font>"), 
            choices = c("no reference population",
                          "1 disease free population",
                          "1 infected population",
                          "1 infected and 1 disease free populations")),
hr(),
HTML("Number of parameters to be evaluated:"),
textOutput("para", inline=TRUE),
p(),
HTML("Degree of freedom:"),
textOutput("ddl", inline=TRUE),
uiOutput("UI"),
hr(),

HTML("<font color='#000080'><b>6: Click the button!</b></font> <br>"),

conditionalPanel(condition = 'input.qtest != 0', 
                 actionButton("button","Estimate")
                ),
hr(),

uiOutput("sortie1"),

hr()
  )
) #End shinyUI
      