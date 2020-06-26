## ----echo=FALSE,results='hide',message=FALSE----------------------------------
######################################################################
### Author: Michael Höhle <http://www.math.su.se/~hoehle>
### Course: The Mathematics and Statistics of Infectious Disease Outbreaks
###         (MT3002 - Summer 2020) at the Department of Mathematics,
###         Stockholm University.
###         (https://kurser.math.su.se/course/view.php?id=911)
###         given by Tom Britton and Michael Höhle
###    
### License:
### This work is licensed under a <a rel="license"
### href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
### Attribution-ShareAlike 4.0 International License</a>.
###
### Description:
###  Rnw File for lecture 05
###
### History:
###  -- 22 Jun 2020 file created 
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(RColorBrewer)
library(lubridate)


## ----fig.show="hide", results="hide"------------------------------------------
# Confirm theory by simulation
beta <- 2.5
gamma <- 1
D <- rexp(1e6, gamma)
C <- rpois(1e6, D*beta)
hist(C, nclass=100)
p <- 1/(beta/gamma +1)
points(0:80, 1e6*dgeom(0:80, prob=p))


## -----------------------------------------------------------------------------
# Parameters
Rt <- 2.5
k  <- 0.45 

# Evaluate on a larger enough grid, so E(Y_t) is determined accurate enough
# We also include -1 in the grid to get a point (0,0) needed for the Lorenz curve
df <- data.frame(x=-1:250) %>% 
  mutate(pmf_negbin= dnbinom(x, mu=Rt, size=k)) %>% 
  mutate(pmf_geom= dgeom(x, p=p))
ggplot(df, aes(x=x, ymin=0, ymax=pmf_negbin, color="negbin")) + 
  geom_linerange() + 
  geom_linerange(aes(x=x + 0.4, ymin=0, ymax=pmf_geom, color="geometric")) +
  scale_color_brewer(type="qual",palette="Set1", name="Offspring:") +
  xlab("Number of secondary cases") + ylab("PMF") + coord_cartesian(xlim=c(0,20)) +
  ggtitle(str_c("R(t)=", Rt, ". k=", k))


## ----simoutbreak, warning=FALSE-----------------------------------------------
####################################################################
#' Simulate time series with this generation time/serial interval and R_e(t)
#' @param n Number of time periods of the outbreak
#' @param Ret A function Re(t) which given t returns the current effective reproduc number
#' @param GT_obj
#' @param Initial number of cases for t=1
#' @return A names vector with the (expected) number of new symptom onsets per day. The outbreak is aligned such that $t=1$ corresponds to 2020-02-15.
####################################################################

routbreak <- function(n=100, Ret, GT_obj, initial_cases = 1) {
  # Set up time series of incident cases
  y <- rep(0, n + length(GT_pmf))
  y[seq_len(length(initial_cases))] <- initial_cases
  # Outbreak starts on 2020-02-15
  dates <- as.Date("2020-02-15") + 0:(n-1)
  # Extract serial interval PMF, ignore support at 0.
  GT_pmf <- GT_obj$GT[-1]
  
  # Loop over all time points
  for (i in 1:n) {
    date <- dates[i]
    y[i + 1:length(GT_pmf)] <- y[i] * (Ret(date) * GT_pmf) + y[i + 1:length(GT_pmf)]
  }
  
  # Data frame with the result. Assume we start on 15th of Feb
  res <- data.frame(Date=dates, y=y[1:n])
  
  #Done
  return(res)
}

# Define time varying effective reproduction number
Ret <- function(date) ifelse(date <= as.Date("2020-03-15"), 2.5, 0.95)
# Define generation time to use
GT_pmf <- structure( c(0, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1), names=0:7)
GT_obj <- R0::generation.time("empirical", val=GT_pmf)
GT_obj

# Generate an outbreak (no stochasticity, just the difference equation)
out <- routbreak(n=60, Ret=Ret, GT_obj=GT_obj)
out <- out %>% mutate(ratio = y/lag(y))
# Data frame with the true values
Ret_true <- data.frame(Date=out$Date) %>% mutate(R_hat=Ret(Date), Method="Truth")

## ----plotsimoutbreak, warning=FALSE-------------------------------------------
p1 <- ggplot(out, aes(x=Date, y=y)) + geom_line() + ylab("No. cases")
p2 <- ggplot(out, aes(x=Date, y=ratio)) + geom_line() + ylab("Ratio")
gridExtra::grid.arrange(p1, p2, ncol=2)


## ----epiestim, echo=TRUE, warning=FALSE---------------------------------------
library(EpiEstim)
# Rename data.frame columns to names handled by the EpiEstim pkg.
out_epiestim <- out %>% rename(I = y, dates = Date) %>% select(dates, I)

# Estimate the instantaneous reproduction number
res <- EpiEstim::estimate_R(out_epiestim, method = "non_parametric_si",
                            config=make_config(si_distr=GT_obj$GT,
                                               t_start=2:nrow(out_epiestim),
                                               t_end=2:nrow(out_epiestim))
)

# Convert result to a data.frame
rt_irt_df <- data.frame(Date=res$dates[res$R$t_end], 
                        R_hat=res$R$`Mean(R)`, 
                        lower=res$R$`Quantile.0.025`, 
                        upper=res$R$`Quantile.0.975`, 
                        Method="R(t)")


## ----PLOTINSTANTANEOUSR0, warning=FALSE, message=FALSE------------------------
ggplot(rt_irt_df, aes(x=Date, y=R_hat, color=Method)) +  
  geom_ribbon(aes(x=Date,  ymin=lower, ymax=upper, color=NULL), alpha=0.1) +
  geom_line(data=Ret_true, aes(x=Date, y=R_hat, color=Method), lwd=2) + 
  geom_line(lwd=1) +
  coord_cartesian(ylim=c(0, 5)) +
  scale_color_brewer(type="qual",palette="Set1", name="Method:") +
  ylab(expression(R(t))) +
  theme_minimal() + 
  theme(legend.position="bottom")


## -----------------------------------------------------------------------------
# Estimate effective reproduction number R(t) for Sweden
# For a detailed report see, e.g., https://www.folkhalsomyndigheten.se/contentassets/4b4dd8c7e15d48d2be744248794d1438/riket-skattning-av-effektiva-reproduktionsnumret-2020-06-12.pdf

library(tidyverse)
library(readxl)

# Download data - note the description at https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/bekraftade-fall-i-sverige/
# which includes the following part:
# 
#  Cases: Antalet fall som redovisas i statistiken är baserat på laboratoriebekräftade fall anmälda enligt smittskyddslagen och redovisas enligt rapporteringsdatum, inklusive positiva prover tagna inom sentinelprovtagningen. 
#   ....
#  Deaths: Tidserien med antalet avlidna per dag innehåller endast de fall där datum för dödsfallet är känt, därför kan totalt antal avlidna skilja sig mot det antal som redovisas i tidsserien. Socialstyrelsen följer antalet avlidna med angiven dödsorsak covid-19. Läs mer om olika datakällor för avlidna i covid-19 i Socialstyrelsen faktablad.
file_name <- str_c("Folkhalsomyndigheten-COVID19-", Sys.Date(), ".xlsx")
if (!file.exists(file_name)) {
  download.file("https://www.arcgis.com/sharing/rest/content/items/b5e7488e117749c19881cce45db13f7e/data", destfile=file_name)
}
covid19_reports <- readxl::read_xlsx(file_name, sheet = "Antal per dag region")
covid19_icu     <- readxl::read_xlsx(file_name, sheet = "Antal intensivvårdade per dag")
# Import of death dates appears not to work, because of "Uppgift saknas in the row" first date is 2020-03-11
# so we do it in two sweeps. 1) read data to determine number of rows 2) read again ignoring the last row. 
# Warning: This can probably be done more elegantly...
covid19_deaths  <- readxl::read_xlsx(file_name, sheet = "Antal avlidna per dag")
covid19_deaths  <- readxl::read_xlsx(file_name, sheet = "Antal avlidna per dag", range=cell_rows(c(1, nrow(covid19_deaths)-1)))

# I'll use values from
# https://www.drugs.com/medical-answers/covid-19-symptoms-progress-death-3536264/
# to define the lags
ts_df <- left_join(covid19_reports %>% select(Statistikdatum, Totalt_antal_fall), 
                   covid19_icu, by=c("Statistikdatum"="Datum_vårdstart")) %>% 
    left_join(covid19_deaths, by=c("Statistikdatum"="Datum_avliden")) %>% 
  rename(Date = Statistikdatum) %>% 
  mutate(ratio_icu_case = Antal_intensivvårdade / lag(Totalt_antal_fall, n=12)) %>% #12
  mutate(ratio_death_icu = Antal_avlidna / lag(Antal_intensivvårdade, n=7)) %>% #7
  mutate(ratio_death_case = Antal_avlidna / lag(Totalt_antal_fall, n=19)) %>% #19
  mutate(Date = as.Date(Date, tc="UTC"))

# Helper function to calculate a centered running mean
roll7 <- function(x) { zoo::rollmean(x, k=7, align="center", na.pad=TRUE) }

# Compute running means for various quantities, replace this code by mutata_at
ts_df <- ts_df %>%
    mutate(Totalt_antal_fall7 = roll7(Totalt_antal_fall),
           Antal_intensivvårdade7 = roll7(Antal_intensivvårdade),
           Antal_avlidna7 = roll7(Antal_avlidna),
           ratio_icu_case7 = roll7(ratio_icu_case),
           ratio_death_icu7 = roll7(ratio_death_icu),
           ratio_death_case7 = roll7(ratio_death_case)
           )

## ----results="asis"-----------------------------------------------------------
# Show some summaries
ts_df %>% 
  select(Totalt_antal_fall, Antal_intensivvårdade, Antal_avlidna) %>% 
  summarise_all(sum, na.rm=TRUE) %>% 
  xtable::xtable(digits=0) %>% 
  print(include.rownames=FALSE)


## ----TS, warning=FALSE--------------------------------------------------------
## Show resulting time series
p1 <- ggplot(ts_df, aes(x=Date, y=Totalt_antal_fall)) + 
  geom_col() +
  geom_line(aes(y=Totalt_antal_fall7), lwd=1.2)
p2 <- ggplot(ts_df, aes(x=Date, y=Antal_intensivvårdade)) + 
  geom_col() + 
  geom_line(aes(y=Antal_intensivvårdade7), lwd=1.2)
p3 <-ggplot(ts_df, aes(x=Date, y=Antal_avlidna)) + 
  geom_col() +
  geom_line(aes(y=Antal_avlidna7), lwd=1.2)

## Show the three plots
gridExtra::grid.arrange(p1, p2, p3, ncol=3)


## ----warning=FALSE------------------------------------------------------------
# Plot of 
p1 <- ggplot(ts_df, aes(x=Date, y=ratio_icu_case)) + geom_line() +  
  geom_line(aes(y=ratio_icu_case7), col="steelblue", lwd=1.2) +
  coord_cartesian(xlim=c(as.Date("2020-04-01"), NA), ylim=c(0,0.3)) +
  ggtitle("ICU(t) / Cases(t-12)")
p2 <- ggplot(ts_df, aes(x=Date, y=ratio_death_icu)) + geom_line() + 
  geom_line(aes(y=ratio_death_icu7), col="steelblue", lwd=1.2) +
  coord_cartesian(xlim=c(as.Date("2020-04-01"), NA), ylim=c(0,5)) +
  ggtitle("Death(t) / ICU(t-7)")
p3 <- ggplot(ts_df, aes(x=Date, y=ratio_death_case)) + geom_line() + 
  geom_line(aes(y=ratio_death_case7), col="steelblue", lwd=1.2) +
  coord_cartesian(xlim=c(as.Date("2020-04-01"), NA), ylim=c(0,1)) +
  ggtitle("Death(t) / Cases(t-19)")

gridExtra::grid.arrange(p1, p2, p3, ncol=3)


## ----warning=FALSE------------------------------------------------------------
#########################################
## R(t) estimation
#########################################

library(EpiEstim)

# Rename data.frame columns to names handled by the EpiEstim pkg.
out_epiestim <- ts_df %>% mutate(Date= as.Date(Date)) %>% 
  rename(I = Totalt_antal_fall, dates = Date) %>% select(dates, I) %>% 
  # Ensure all reports are there (6 days back, this is rather conservative)
  filter(dates <= Sys.Date()-6) 

# From the FOHM report: Här har vi använt det skattade serieintervallet från Nishiura et al.
# (2020) med ett medelvärde på 4.8 dagar och en standaravvikelse på 2.3 dagar.
# Note: Weibull distribution in the Nishiura et al. paper vs. (shifted) gamma in EpiEstim

# Estimate the instantaneous R_0 - no smoothing
res1 <- EpiEstim::estimate_R(out_epiestim, method = "parametric_si",
                            config=make_config(mean_si = 4.8, std_si = 2.3,
                                               t_start=2:nrow(out_epiestim),
                                               t_end=2:nrow(out_epiestim)))
# Use smoothing over 7 days (recommended!)
res7 <- EpiEstim::estimate_R(out_epiestim, method = "parametric_si",
                            config=make_config(mean_si = 4.8, std_si = 2.3,
                                               t_start=2:(nrow(out_epiestim)-6),
                                               t_end=8:nrow(out_epiestim)))

# Convert result to a data.frame
epiestim2df <- function(res) {
  data.frame(Date=res$dates[res$R$t_end], R_hat=res$R$`Mean(R)`, lower=res$R$`Quantile.0.025(R)`, upper=res$R$`Quantile.0.975`)
}

rt_irt7_df <- epiestim2df(res7) %>% mutate(method = "smooth7")
rt_irt1_df <- epiestim2df(res1) %>% mutate(method = "no smoothing")
rt_df <- rbind(rt_irt7_df, rt_irt1_df)

## ----rtplots, warning=FALSE---------------------------------------------------
# The plot - doesn't make too much sense because testing is substantially 
# increased in the last weeks
# Start at 2020-03-16 so imported casesd do not play a role anymore
ggplot(rt_df, aes(x=Date, y=R_hat)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="steelblue", alpha=0.2) +
  geom_line() +
  scale_x_date(breaks = seq(as.Date("2020-03-16"), as.Date("2020-06-30"),by="1 week"), date_labels = "%b %d") +
  coord_cartesian(xlim=c(as.Date("2020-03-20"), NA), ylim=c(0, 2)) +
  ylab(expression(R(t)))  +
  geom_hline(yintercept=1, lty=2, col="salmon2") +
  facet_wrap(~ method, ncol=1)


## ----scrape_tests-------------------------------------------------------------
#######
# Scrape the weekly FOHM data on testing in order to compute proportion 
#######

library(rvest)
# Extract table with the number of test results results from round 1
results <- read_html("https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/antal-individer-som-har-testats-for-covid-19/")

tests <- results %>%
  html_nodes(css = "table , caption") %>%
  .[[1]] %>%
  html_table(header = 1) %>%
  as_tibble()

# As table
tests <- tests %>%
  head(n = -1) %>%
  mutate(week = str_replace(Vecka, "vecka ", ""), 
         n_tested = str_replace(`Testade individer`," ", "")) %>% mutate_at(c("week", "n_tested"), as.numeric) %>% select(week, n_tested)

ts_week <- ts_df %>% mutate(week = week(Date)) %>% group_by(week) %>% 
  summarise(n_positive = sum(Totalt_antal_fall, na.rm=TRUE))

tab <- inner_join(tests, ts_week, by=c("week")) %>% mutate(prop_positive = n_positive / n_tested)

## ----plottests, fig.height=4--------------------------------------------------
p1 <- ggplot(tab, aes(x = as.integer(week), y = n_tested)) +
  geom_col() +
  ylim(0, NA) +
  xlab("Week") +
  scale_x_continuous(breaks = min(tab$week):max(tab$week)) +
  theme_minimal()
p2 <- ggplot(tab, aes(x = as.integer(week), y = n_positive / n_tested)) +
  geom_line() +
  xlab("Week") +
  scale_y_continuous(labels = scales::percent, limit = c(0, NA)) +
  scale_x_continuous(breaks = min(tab$week):max(tab$week)) +
  theme_minimal()
gridExtra::grid.arrange(p1,p2, ncol=2)


## ----calcrtgermany------------------------------------------------------------
# Load newest data from the RKI web-page, see https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/R-Wert-Erlaeuterung.pdf?__blob=publicationFile for details, which
# is part of https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/Nowcasting.html
file_name <- str_c("Nowcasting_Zahlen-",Sys.Date(),  ".xlsx")
if (!file.exists(file_name)) {   
  file_url <- "https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/Nowcasting_Zahlen.xlsx?__blob=publicationFile"
  download.file(url=file_url,destfile= file_name, mode="wb")
}
# Load Excel-File
data <- xlsx::read.xlsx(file = file_name, sheetName = "Nowcast_R",  encoding = "UTF-8" ) 

# Rename data.frame columns to names handled by the EpiEstim pkg.
out_epiestim <- data %>% 
  rename(dates=Datum.des.Erkrankungsbeginns, 
         I=Punktschätzer.der.Anzahl.Neuerkrankungen..ohne.Glättung.)  %>% 
  select(dates, I)

# Use smoothing over 7 days (recommended!)
res7 <- EpiEstim::estimate_R(out_epiestim, method = "parametric_si",
                            config=make_config(mean_si = 4.8, std_si = 2.3,
                                               t_start=2:(nrow(out_epiestim)-6),
                                               t_end=8:nrow(out_epiestim)))

rt_df <- epiestim2df(res7) %>% mutate(method = "smooth7")

## ----rtgermany----------------------------------------------------------------
p1 <- ggplot(out_epiestim, aes(x=dates, y=I)) + geom_col() + 
  xlab("Disease onset (estimated)") + ylab("Counts") +   
  scale_x_date(breaks = seq(as.Date("2020-03-02"), as.Date("2020-06-30"),by="1 week"), date_labels = "%b %d") +
  coord_cartesian(xlim=c(as.Date("2020-03-02"), NA)) 

p2 <- ggplot(rt_df, aes(x=Date, y=R_hat)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="steelblue", alpha=0.2) +
  geom_line() +
  scale_x_date(breaks = seq(as.Date("2020-03-02"), as.Date("2020-06-30"),by="1 week"), date_labels = "%b %d") +
  coord_cartesian(xlim=c(as.Date("2020-03-02"), NA), ylim=c(0, 2)) +
  ylab(expression(R(t)))  +
  geom_hline(yintercept=1, lty=2, col="salmon2") +
  facet_wrap(~ method, ncol=1)

gridExtra::grid.arrange(p1,p2, ncol=1)

