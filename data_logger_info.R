Sys.setlocale("LC_TIME", "C")

temp_data <- 
  read.table('labeglo2/data_from_logger/ibutton_lake14_20192020.tsv', 
             header = F)
temp_data_summer <- read.table('labeglo2/data_from_logger/2019_lake14.csv', 
                               header = T, sep = ',')
colnames(temp_data_summer) <- c('day', 'temperature')
temp_data_summer$day <- sub(' .*', '', temp_data_summer$day)

colnames(temp_data) <- c('day', 'time', 'temperature')
temp_data$temperature <- sub(',', '.', temp_data$temperature)
temp_data$temperature <- as.numeric(temp_data$temperature)
temp_data$day <- sub('\\.20', '\\.', temp_data$day)

temp_data <- temp_data[-2]
temp_data_year <- rbind(temp_data_summer, temp_data)
temp_data_day_average <- aggregate(temperature ~ day, data = temp_data_year, mean)

temp_data_day_average$day <- 
  as.POSIXct(as.character(temp_data_day_average$day), '%d.%m.%y', tz = " ")
temp_data_day_average <- 
  temp_data_day_average[temp_data_day_average$day!=min(temp_data_day_average$day),]

date <- as.Date(temp_data_day_average$day)
d <- data.frame(date, temp_data_day_average$temperature, 
                month = format(date, '%B'))

ggplot(temp_data_day_average, aes(day, temperature)) +
  geom_vline(xintercept = as.POSIXct("2019-09-19", '%y-%m-%d'), color = 'grey50') +
  geom_vline(xintercept = as.POSIXct("2019-11-06", '%y-%m-%d'), color = 'grey50') +
  geom_vline(xintercept = as.POSIXct("2019-12-09", '%y-%m-%d'), color = 'grey50') +
  geom_vline(xintercept = as.POSIXct("2020-01-25", '%y-%m-%d'), color = 'grey50') +
  geom_vline(xintercept = as.POSIXct("2020-06-20", '%y-%m-%d'), color = 'grey50') +
  geom_hline(yintercept = 1.5, color = 'firebrick', alpha = 0.2) +
  geom_point(size = 0.2) +
  geom_line() +
  scale_y_continuous("Mean day temperature", breaks = c(0:12), 
                                             limits = c(0,12)) +
  xlab('Date') +
  #scale_x_date(name = "Date", 
               #date_labels = "%e", 
  #             breaks = as.Date(temp_data_day_average$day,"%y-%m-%d"),
  #             date_breaks = "5 days") +
  #scale_x_date(name = "Date", breaks = as.Date(temp_data_day_average$day, 
  #                                             "%y-%m-%d")) +
  theme_light()

ggsave('labeglo2/data_from_logger/lake14_2019-2020_1m.png')  
  
  
  