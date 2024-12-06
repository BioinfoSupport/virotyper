library(bslib)
library(shiny)
library(bsicons)
library(tidyverse)
future::plan(future::multisession)

# Accept up to 5Gb file upload
options(shiny.maxRequestSize=5000*1024^2)

TerminalTask <- R6::R6Class(
  "TerminalTask",
  inherit = shiny::ExtendedTask,
  public = list(
    fifo = NULL, # FIFO where stdin and stdout are redirected
    fifo_con = NULL,
    lines = character(0), # Terminal content (what has been read from the fifo)
    
    initialize = function() {
      self$fifo <- tempfile(fileext = ".fifo") 
      self$fifo_con <- fifo(self$fifo,"w+") # Create the FIFO
      message(self$fifo)
      super$initialize(function(cmd) { # Call superclass
        cmd <- stringr::str_glue("({cmd}) > {self$fifo} 2>&1") # Adapt command line to redirect output to the FIFO
        promises::future_promise(system(cmd,intern=FALSE))
      })
    },
    
    finalize = function() {
      close(self$fifo_con)
    },
    
    # Get current terminal content
    current_content = function() {
      self$lines <- c(self$lines,readLines(self$fifo_con))
      self$lines
    }
  )
)


init_workdir <- function(workdir=tempfile()) {
  message(workdir)
  
  # Copy data/db
  dir.create(file.path(workdir,"data/db"),recursive = TRUE)
  file.copy("data/db",file.path(workdir,"data"),recursive = TRUE)

  # Copy bin/
  file.copy("bin",workdir,recursive = TRUE)

  # Copy notebooks/
  file.copy("notebooks",workdir,recursive = TRUE)
  
  # Copy Makefile
  file.copy("Makefile",workdir)

  #print(list.files(workdir,recursive = TRUE,full.names = TRUE))
  workdir
}

ui_panel_db <- function() {
  div(class="shiny-input-panel",
    selectInput("resistance_db",list("Select Resistance DB (or ",actionLink("show_upload_db_dialog","upload a new DB"),")"),choices=character(0),width="100%")
  )
}


ui_panel_vcf <- function() {
  nav_panel(
    title = "From VCF",
    tags$p("Use this form to generate a resistance report from a VCF file containing already called mutation against the reference database."),
    div(class="shiny-input-panel",
      fileInput("vcf_file","Upload a VCF file",accept = c(".vcf",".gz"),multiple = FALSE,placeholder = "Upload a VCF file (with required INFO fields DP,AF,BCSQ)",width = "100%")
    ),
    input_task_button("btn_vcf_cmd_execute","Run"),
    shinyjs::disabled(actionButton("btn_view_vcf_report","View HTML report",icon = icon("eye"))),
    shinyjs::disabled(downloadButton("btn_dl_vcf_report_docx","Download DOCX Report",icon = icon("download"))),
    br(),
    verbatimTextOutput("vcf_term_output")
  )
}

ui_panel_fasta <- function() {
  nav_panel(
    title = "From FASTA",
    tags$p("Use this form to generate a resistance report from a FASTA file containing an assembly. The assembly will be mapped to against the reference database, and variant calling performed."),
    div(class="shiny-input-panel",
        fileInput("fasta_file","Upload a FASTA file",accept = ".fasta",multiple = FALSE,placeholder = "Upload a FASTA file containing your assembly",width = "100%")
    ),
    verbatimTextOutput("fasta_term_output")
  )
}

ui <- function() {
  page_fluid(
    title = ,
    theme = bs_theme(),
    titlePanel(list(icon("biohazard"),"Virotyper"),"Virotyper"),
    shinyjs::useShinyjs(),

    ui_panel_db(),
    br(),
    navset_tab(
      nav_panel(
        "Viral Resistance Analysis",
        p("Generate a viral resistance report"),
        navset_tab(
          ui_panel_vcf(),
          ui_panel_fasta()
        )
      )
    )
    
    #   selectInput("resistance_db","Resistance database",choices = c("HSV1","HSV2")),
    #   textAreaInput("fasta_seq","DNA sequences (FASTA)",placeholder = ">UL23\nATGCTCATGACTGCATCTA...\n>UL30\nATGCTCATGACTGCATCTA...",height = "400px"),
    #   shinycssloaders::withSpinner(textOutput("terminal_stdout",container = pre),type = 1)
  )
}

server <- function(input, output, session) {

  # Initialize a new working directory for the current session
  workdir <- init_workdir()

  # Startup Disclaimer Message
  showModal(modalDialog(easyClose = FALSE,title = "Disclaimer",footer = modalButton("Accept"),class="bg-warning",
     card_body("IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
        ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
        TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
        OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
     )
  ))
    
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Resistance Database Selection and Upload logic
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  available_dbs <- reactiveVal(list.files(file.path(workdir,"data/db")))

  observeEvent(available_dbs(),{
    updateSelectInput(inputId="resistance_db",choices=available_dbs())
  })

  # When upload DB link is clicked: show the modal box
  observeEvent(input$show_upload_db_dialog,{
    showModal(modalDialog(
      title = "DB upload",
      fileInput("db_file","TAR achive of the new DB",accept = ".gz",multiple = FALSE,placeholder = "Upload a resistance database archive (.tar.gz)"),
      easyClose = TRUE,size = "xl",
      footer = NULL
    ))
  })
  
  observeEvent(input$db_file,{
    validate(need(is.character(input$db_file$datapath),"invalid TAR"))
    message("untaring input$db_file$datapath")
    withProgress(message="Uncompress TAR database",{
      untar(input$db_file$datapath,exdir=file.path(workdir,"data/db"))
    })
    available_dbs(list.files(file.path(workdir,"data/db")))
    removeModal()
  })
  
  
  
  #-#-#-#-#-#-#-#-#-#-#
  # VCF panel logic
  #-#-#-#-#-#-#-#-#-#-#
  vcf_task <- TerminalTask$new() |>
    bind_task_button("btn_vcf_cmd_execute")
  
  current_vcf <- reactiveVal(NULL)
  
  # When a new VCF file is uploaded
  observeEvent(input$vcf_file,{
    validate(need(is.character(input$vcf_file$datapath),"invalid VCF"))
    # Copy inputs to workdir location
    dir.create(file.path(workdir,"input"),recursive = TRUE)
    current_vcf <- file.path("input",input$vcf_file$name)
    file.copy(input$vcf_file$datapath, file.path(workdir,current_vcf), overwrite = TRUE)
    current_vcf(current_vcf)
  })
  
  # When Run is clicked
  observeEvent(input$btn_vcf_cmd_execute,{
    # Generate report
    cmd <- str_glue("echo {workdir} && cd '{workdir}' && make DB_ID={input$resistance_db} '{current_vcf()}.all'")
    message(cmd)
    vcf_task$invoke(cmd)
  })
  
  # When View Button is clicked
  observeEvent(input$btn_view_vcf_report,{
    showModal(modalDialog(
      title = list("VCF report",downloadButton("btn_dl_vcf_report_html","Download",icon = icon("download"),class="btn-primary")),
      tags$iframe(srcdoc=readLines(str_glue("{file.path(workdir,current_vcf())}.{input$resistance_db}.html")), style='width:100%;height:600px;'),
      easyClose = TRUE,size = "xl",footer=list()
    ))
  })
  
  # When running task status change
  observe({
    switch (vcf_task$status(),
      "success" = {
        shinyjs::enable("btn_view_vcf_report")
        shinyjs::enable("btn_dl_vcf_report_docx")
        shinyjs::enable("vcf_file")
      },
      "running" = {
        shinyjs::disable("btn_view_vcf_report")
        shinyjs::disable("btn_dl_vcf_report_docx")
        shinyjs::disable("vcf_file")
      }
    )
  })
  
  # When input parameters change
  observeEvent(input$vcf_file,{
    shinyjs::disable("btn_view_vcf_report")
    shinyjs::disable("btn_dl_vcf_report_docx")
  })
  
  # When the download HTML button is clicked 
  output$btn_dl_vcf_report_html <- downloadHandler(
    filename = function(){str_glue("{basename(file.path(workdir,current_vcf()))}.{input$resistance_db}.html")},
    content = function(file) {file.copy(str_glue("{file.path(workdir,current_vcf())}.{input$resistance_db}.html"),file)}
  )
  
  # When the download DOCX button is clicked
  output$btn_dl_vcf_report_docx <- downloadHandler(
    filename = function(){str_glue("{basename(file.path(workdir,current_vcf()))}.{input$resistance_db}.docx")},
    content = function(file) {file.copy(str_glue("{file.path(workdir,current_vcf())}.{input$resistance_db}.docx"),file)}
  )
  
  # Periodic update of Terminal output
  output$vcf_term_output <- renderText({
    if (vcf_task$status() == "running") invalidateLater(500)
    paste0(vcf_task$current_content(),collapse = "\n")
  })

}

shinyApp(ui = ui(),server = server)
