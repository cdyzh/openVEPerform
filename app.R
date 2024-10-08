library(shiny)
library(dplyr)
library(yogiroc)
library(httr)
library(jsonlite)
library(shinycssloaders)
library(DT)

# Define UI for the application
ui <- fluidPage(
  titlePanel("openVEPerform"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("input_type", "Select Input Type:",
                   choices = list("Chromosome, Position, Ref, Alt" = "chrom_pos",
                                  "Transcript ID and HGVSC" = "hgvsc")),
      
      # File input for Chromosome, Position, Reference_Base, Alternate_Base
      conditionalPanel(
        condition = "input.input_type == 'chrom_pos'",
        fileInput("variant_file", "Upload CSV File (Chromosome, Position, Reference, Alternate)",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
      ),
      
      # File input for Transcript ID and HGVSC
      conditionalPanel(
        condition = "input.input_type == 'hgvsc'",
        fileInput("variant_file_hgvsc", "Upload CSV File (Transcript ID, HGVSC)",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
      ),
      
      # Option to run OpenCRAVAT and fetch the scores
      actionButton("fetchButton", "Fetch VEP Data"),
      
      # Additional controls for the plot
      checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
      checkboxGroupInput("scores", "Select Scores to Include:",
                         choices = list(
                           "VARITY",
                           "REVEL",
                           "AlphaMissense"
                         ),
                         selected = c("VARITY", "REVEL", "AlphaMissense")
      ),
      actionButton("plotButton", "Generate PRC Plot"),
      downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
      downloadButton("downloadPlotPDF", "Download PRC Plot and Metadata"),
      helpText("Upload a CSV file with columns: Chromosome, Position, Reference_Base, Alternate_Base, or Transcript ID and HGVSC.")
    ),
    
    mainPanel(
      plotOutput("prcPlot", width = "600px", height = "600px"),
      textOutput("errorText")
    )
  )
)

server <- function(input, output, session) {
  variant_data <- reactiveVal(NULL)
  plot_info <- reactiveVal(NULL)
  
  observeEvent(input$fetchButton, {
    req(input$input_type)
    
    # Initialize empty list to store results
    all_results <- list()
    
    if (input$input_type == "chrom_pos") {
      req(input$variant_file)
      df <- read.csv(input$variant_file$datapath, stringsAsFactors = FALSE)
      colnames(df) <- trimws(colnames(df))
      
      # Check for necessary columns
      if (!all(c("Chromosome", "Position", "Reference_Base", "Alternate_Base") %in% colnames(df))) {
        output$errorText <- renderText("Error: The uploaded CSV must contain the columns 'Chromosome', 'Position', 'Reference_Base', 'Alternate_Base'.")
        return(NULL)
      }
    } else if (input$input_type == "hgvsc") {
      req(input$variant_file_hgvsc)
      df_hgvsc <- read.csv(input$variant_file_hgvsc$datapath, stringsAsFactors = FALSE)
      colnames(df_hgvsc) <- trimws(colnames(df_hgvsc))
      
      # Ensure the necessary columns exist and rename to standard column name
      if (!all(c("Transcript_ID", "HGVSC") %in% colnames(df_hgvsc))) {
        output$errorText <- renderText("Error: The uploaded CSV must contain the columns 'Transcript_ID' and 'HGVSC'.")
        return(NULL)
      }
      
      # Concatenate the Transcript ID and HGVSC into a single column
      df_hgvsc <- df_hgvsc %>%
        mutate(OpenCRAVAT_Input = paste0(Transcript_ID, ":", HGVSC))
      
      # Use the concatenated column for further processing
      df <- df_hgvsc %>% select(OpenCRAVAT_Input)
    }
    
    # Define the standard column names to ensure consistency
    standard_colnames <- c(
      "base__hugo", 
      "base__achange", 
      "gnomad__af", 
      "varity_r__varity_r", 
      "revel__score", 
      "alphamissense__am_pathogenicity", 
      "classification"
    )
    
    # Show a progress bar while fetching data
    withProgress(message = 'Fetching Variant Data', value = 0, {
      for (i in 1:nrow(df)) {
        if (input$input_type == "chrom_pos") {
          # Check if 'Chromosome' already contains the 'chr' prefix
          if (grepl("^chr", df$Chromosome[i])) {
            chrom <- df$Chromosome[i]  # Use as is if it contains 'chr'
          } else {
            chrom <- paste0("chr", df$Chromosome[i])  # Add 'chr' prefix if missing
          }
          pos <- df$Position[i]
          ref <- df$Reference_Base[i]
          alt <- df$Alternate_Base[i]
          
          # Construct API URL for chromosome-based input
          api_url <- paste0(
            "https://run.opencravat.org/submit/annotate?",
            "chrom=", chrom,
            "&pos=", pos,
            "&ref_base=", ref,
            "&alt_base=", alt,
            "&annotators=clinvar,gnomad,varity_r,revel,alphamissense"
          )
        } else if (input$input_type == "hgvsc") {
          hgvsc_input <- df$OpenCRAVAT_Input[i]
          
          print(hgvsc_input) # DEBUG
          # Construct API URL for HGVSC input
          api_url <- paste0(
            "https://run.opencravat.org/submit/annotate?",
            "hgvs=", hgvsc_input,
            "&annotators=clinvar,gnomad,varity_r,revel,alphamissense"
          )
        }
        
        # Make the GET request to OpenCRAVAT
        response <- GET(api_url)
        
        # Parse the JSON response
        result <- fromJSON(content(response, "text"), flatten = TRUE)
        
        # Handle cases where certain annotations might be missing
        clinvar_sig <- ifelse(!is.null(result$clinvar$sig), result$clinvar$sig, NA)
        
        # Determine the classification value (T/F) based on ClinVar significance
        classification <- NA
        if (!is.na(clinvar_sig)) {
          clinvar_sig_lower <- tolower(clinvar_sig)
          
          if (grepl("benign", clinvar_sig_lower)) {
            classification <- TRUE
          } else if (grepl("pathogenic", clinvar_sig_lower) && !grepl("conflicting", clinvar_sig_lower)) {
            classification <- FALSE
          }
        }
        
        # Extract into corresponding columns - TODO: be able to add additional predictors
        gene <- ifelse(!is.null(result$crx$hugo), result$crx$hugo, NA)
        achange <- ifelse(!is.null(result$crx$achange), result$crx$achange, NA)
        gnomad_af <- ifelse(!is.null(result$gnomad$af), result$gnomad$af, NA)
        varity_r <- ifelse(!is.null(result$varity_r$varity_r), result$varity_r$varity_r, NA)
        revel_score <- ifelse(!is.null(result$revel$score), result$revel$score, NA)
        alphamissense_path <- ifelse(!is.null(result$alphamissense$am_pathogenicity), result$alphamissense$am_pathogenicity, NA)
        
        # Create a data frame with consistent column names
        result_df <- data.frame(
          base__hugo = gene,
          base__achange = achange,
          gnomad__af = gnomad_af,
          varity_r__varity_r = varity_r,
          revel__score = revel_score,
          alphamissense__am_pathogenicity = alphamissense_path,
          classification = classification,
          stringsAsFactors = FALSE
        )
        
        # Ensure that all columns match the standard names, even if some are missing
        missing_cols <- setdiff(standard_colnames, colnames(result_df))
        result_df[missing_cols] <- NA  # Add missing columns with NA values
        
        # Append to the results list
        all_results[[i]] <- result_df
        
        incProgress(1 / nrow(df))
        showNotification(paste(i, "/", nrow(df), "variants fetched"), duration = 3, type = "message")
      }
    })
    print("done1")
    
    # Combine all results into a data frame
    variant_data_df <- do.call(rbind, all_results)

    # Update reactive variable
    variant_data(variant_data_df)
    print("done2")
    # Clear any error message
    output$errorText <- renderText("")
  })
  
  observeEvent(input$plotButton, {
    df <- variant_data()
    df <- df[!is.na(df$classification), ]
    print(df) # DEBUG
    
    if (is.null(df) || nrow(df) == 0) {
      output$errorText <- renderText("Not enough rows to generate the PRC plot.")
      return()
    }
    
    gene_s <- df$base__hugo[1]
    exclude_common_variants <- input$common_variant_filter
    selected_scores <- input$scores
    
    # Rename columns to match what you need for PRC plotting
    names(df)[names(df) == "varity_r__varity_r"] <- "VARITY"
    names(df)[names(df) == "alphamissense__am_pathogenicity"] <- "AlphaMissense"
    names(df)[names(df) == "revel__score"] <- "REVEL"
    names(df)[names(df) == "gnomad__af"] <- "gnomAD_AF"
    
    if (exclude_common_variants) {
      df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
    }
    
    df <- df[order(df$classification), ]
    
    # Filter for the selected scores
    prcfiltered <- df %>%
      filter(rowSums(!is.na(df[selected_scores])) > 0)
    
    B_org <- sum(prcfiltered$classification == TRUE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
    P_org <- sum(prcfiltered$classification == FALSE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
    
    tryCatch({
      yrobj <- yr2(truth = prcfiltered[["classification"]], scores = prcfiltered[selected_scores], high = rep(FALSE, length(selected_scores)))
      
      # Store plot information in plot_info instead of overwriting variant_data
      plot_info(list(
        yrobj = yrobj,
        lty_styles = c("dashed", "solid", "dashed")[1:length(selected_scores)],
        col_styles = c("purple", "cadetblue2", "orange")[1:length(selected_scores)],
        gene_s = gene_s,
        selected_scores = selected_scores,
        B_org = B_org,
        P_org = P_org,
        prcfiltered = prcfiltered  # Save the filtered data for metadata
      ))
      
      output$errorText <- renderText("")
      
    }, error = function(e) {
      plot_info(NULL)
      output$errorText <- renderText("Not enough data to generate PRC plot.")
    })
    
    output$prcPlot <- renderPlot({
      plot_details <- plot_info()
      if (!is.null(plot_details)) {
        tryCatch({
          draw.prc(plot_details$yrobj, lty = plot_details$lty_styles, col = plot_details$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_details$gene_s, " PRCs for ", paste(plot_details$selected_scores, collapse = ", ")))
          abline(h = 90, lty = "dashed")
          legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_details$P_org), paste("# of Benign and Likely Benign:", plot_details$B_org)), pch = 15, bty = "n")
        }, error = function(e) {
          showModal(modalDialog(
            title = 'Error',
            'Not enough data',
            easyClose = TRUE,
            footer = NULL
          ))
        })
      }
    }, width = 600, height = 600, res = 72)
  })
  
  output$downloadPlotPNG <- downloadHandler(
    filename = function() {
      paste("PRC_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      plot_details <- plot_info()
      if (!is.null(plot_details)) {
        png(file, width = 6, height = 6, units = "in", res = 72)
        draw.prc(plot_details$yrobj, lty = plot_details$lty_styles, col = plot_details$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_details$gene_s, " PRCs for ", paste(plot_details$selected_scores, collapse = ", ")))
        abline(h = 90, lty = "dashed")
        legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_details$P_org), paste("# of Benign and Likely Benign:", plot_details$B_org)), pch = 15, bty = "n")
        dev.off()
      }
    }
  )
  
  output$downloadPlotPDF <- downloadHandler(
    filename = function() {
      paste("PRC_Report_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      plot_details <- plot_info()
      if (!is.null(plot_details)) {
        # Generate the PDF report
        rmarkdown::render(input = "report_template.Rmd",
                          output_file = file,
                          params = list(
                            gene_s = plot_details$gene_s,
                            selected_scores = plot_details$selected_scores,
                            B_org = plot_details$B_org,
                            P_org = plot_details$P_org,
                            prcfiltered = plot_details$prcfiltered
                          ),
                          envir = new.env(parent = globalenv()))
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)