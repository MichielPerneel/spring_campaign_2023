# Load necessary libraries
library(dplyr)
library(stringr)

# Function to extract relevant data from a file
extract_data <- function(file_path) {
  # Read the first 36 lines of the file
  file_content <- suppressWarnings(readLines(file_path, warn = FALSE, n = 36))

  # Extract date and time
  date_line <- grep("Date:", file_content, value = TRUE)
  time_line <- grep("Time:", file_content, value = TRUE)

  # Extract the date and time values
  date_value <- as.Date(str_extract(date_line, "\\w+ \\d{1,2}, \\d{4}"), format = "%b %d, %Y")
  time_value <- str_extract(time_line, "\\d{2}:\\d{2}")

  # Combine date and time into a single datetime value
  datetime_value <- as.POSIXct(paste(date_value, time_value), format = "%Y-%m-%d %H:%M")
  # Add one hour to the datetime to adjust for the time lag
  datetime_value <- datetime_value + 3600  # 3600 seconds = 1 hour

  # Extract the blank value
  blank_line <- grep("Blank:", file_content, value = TRUE)
  blank_value <- as.numeric(str_extract(blank_line, "\\d+"))

  # Extract the rPE fit values
  alpha_line <- grep("^Alpha:", file_content, value = TRUE)
  beta_line <- grep("^Beta:", file_content, value = TRUE)
  ek_line <- grep("^Ek:", file_content, value = TRUE)
  ekbeta_line <- grep("^EkBeta:", file_content, value = TRUE)
  rpm_line <- grep("^rPm:", file_content, value = TRUE)
  jvpIIm_line <- grep("^JVPIIm:", file_content, value = TRUE)
  gopIIm_line <- grep("^GOPIIm:", file_content, value = TRUE)

  # Extract numeric values from the lines
  alpha_value <- as.numeric(str_extract(alpha_line, "-?\\d+\\.?\\d*"))
  beta_value <- as.numeric(str_extract(beta_line, "-?\\d+\\.?\\d*"))
  ek_value <- as.numeric(str_extract(ek_line, "-?\\d+\\.?\\d*"))
  ekbeta_value <- as.numeric(str_extract(ekbeta_line, "-?\\d+\\.?\\d*"))
  rpm_value <- as.numeric(str_extract(rpm_line, "-?\\d+\\.?\\d*"))
  jvpIIm_value <- as.numeric(str_extract(jvpIIm_line, "-?\\d+\\.?\\d*"))
  gopIIm_value <- as.numeric(str_extract(gopIIm_line, "-?\\d+\\.?\\d*"))

  # Replace negative values with zero
  alpha_value <- ifelse(alpha_value < 0, 0, alpha_value)
  beta_value <- ifelse(beta_value < 0, 0, beta_value)
  ek_value <- ifelse(ek_value < 0, 0, ek_value)
  ekbeta_value <- ifelse(ekbeta_value < 0, 0, ekbeta_value)
  rpm_value <- ifelse(rpm_value < 0, 0, rpm_value)
  jvpIIm_value <- ifelse(jvpIIm_value < 0, 0, jvpIIm_value)
  gopIIm_value <- ifelse(gopIIm_value < 0, 0, gopIIm_value)

  # Estimate Primary production by multiplying the carbon fixation rate (GOP)
  # with the molar mass of carbon (12 g/mol)
  PP_value <- gopIIm_value * 12

  # Extract the filename
  filename <- basename(file_path)

  # Extract the station number from the filename
  station <- as.numeric(str_extract(filename, "(?<=Spring_)[0-9]+"))
  # Extract the sample number from the filename (after 'st')
  sample <- as.numeric(str_extract(filename, "(?<=st)\\d+"))

  # Combine all extracted data into a dataframe
  data_frame <- data.frame(
    Filename = filename,
    Station = station,
    Sample = sample,
    Datetime = datetime_value,
    Blank = blank_value,
    Alpha = alpha_value,
    Beta = beta_value,
    Ek = ek_value,
    EkBeta = ekbeta_value,
    rPm = rpm_value,
    JVPIIm = jvpIIm_value,
    GOPIIm = gopIIm_value,
    PP = PP_value,
    stringsAsFactors = FALSE
  )

  return(data_frame)
}

# Function to process all files in a directory and return a dataframe
process_directory <- function(directory_path) {
  # List all files in the directory
  file_list <- list.files(directory_path, full.names = TRUE, recursive = TRUE)

  # Initialize an empty list to store all results
  all_data_list <- list()

  # Loop through each file and extract data
  for (file in file_list) {
    # Check if the file is a text file
    if (grepl("\\.txt$", file)) {
      # Extract data from the file
      file_data <- extract_data(file)

      # Append the data to the list
      all_data_list <- append(all_data_list, list(file_data))
    }
  }

  # Combine all data frames into one and remove duplicates
  all_data <- do.call(rbind, all_data_list) %>%
    distinct()

  return(all_data)
}

# Directory path for 130
root_directory_130 <- "data/raw/LabSTAF/text_files/130"

# Process the directory and get the dataframe for Station 130
data_130 <- process_directory(root_directory_130)

# Directory path for 51
root_directory_51 <- "data/raw/LabSTAF/text_files/51"

# Process the directory and get the dataframe for Station 51
data_51 <- process_directory(root_directory_51)

# Combine the two datasets into one
final_data <- bind_rows(data_130, data_51)

# Print the combined dataframe for inspection
print(head(final_data))

# Optionally, write the final dataframe to a CSV file
write.csv(final_data, file = "data/raw/LabSTAF/labstaf_combined_data.csv", row.names = FALSE)
