library(bslib)
library(shiny)
library(bsicons)
library(tidyverse)

# Accept up to 5Gb file upload
options(shiny.maxRequestSize=5000*1024^2)


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
  inputPanel(
    title = "Select Resistance Database",
    selectInput("resistance_db",list("Select Resistance DB (or ",actionLink("show_upload_db_dialog","upload a new DB"),")"),choices=character(0),)
  )
}


ui_panel_vcf <- function() {
  nav_panel(
    title = "From VCF",
    tags$p("Use this form to generate a resistance report from a VCF file containing already called mutation against the reference database."),
    inputPanel(
      title = "Generate a report from a VCF containing variants",
      fileInput("vcf_file","VCF file",accept = c(".vcf",".gz"),multiple = FALSE,placeholder = "Upload a VCF file (with INFO fields DP,AF,BCSQ)"),
    ),
    shinyjs::disabled(actionButton("btn_view_vcf_report","View HTML report",icon = icon("eye"))),
    shinyjs::disabled(downloadButton("btn_dl_vcf_report_docx","Download DOCX Report",icon = icon("download"))),
    verbatimTextOutput("term_output")
  )
}

ui <- function() {
  page_fluid(
    title = ,
    theme = bs_theme(),
    titlePanel(list(icon("biohazard"),"Virotyper"),"Virotyper"),
    shinyjs::useShinyjs(),

    ui_panel_db(),    
    navset_tab(
      nav_panel(
        "Viral Resistance Analysis",
        p("Generate a viral resistance report"),
        navset_tab(
          ui_panel_vcf(),
          nav_panel("From FASTA"),
          nav_panel("From BAM")
        )
      ),
      nav_panel(
        "ONT Reads Analysis"
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
  current_vcf <- reactiveVal(NULL)
  
  # When a new VCF file is uploaded
  observeEvent(input$vcf_file,{
    validate(need(is.character(input$vcf_file$datapath),"invalid VCF"))
    shinyjs::disable("btn_view_vcf_report")
    shinyjs::disable("btn_dl_vcf_report_docx")
    
    # Copy inputs to workdir location
    dir.create(file.path(workdir,"input"),recursive = TRUE)
    current_vcf <- file.path(workdir,"input",input$vcf_file$name)
    file.copy(input$vcf_file$datapath, current_vcf, overwrite = TRUE)
    current_vcf(current_vcf)
    
    # Generate report
    cmd <- str_glue("cd '{workdir}' && make DB_ID={input$resistance_db} '{current_vcf}.all'")
    message(cmd)
    withProgress(message="Generate VCF report",{
      system(cmd)
    })
    shinyjs::enable("btn_view_vcf_report")
    shinyjs::enable("btn_dl_vcf_report_docx")
  })
  
  observeEvent(input$btn_view_vcf_report,{
    showModal(modalDialog(
      title = "VCF report",
      tags$iframe(srcdoc=readLines(str_glue("{current_vcf()}.{input$resistance_db}.html")), style='width:100%;height:600px;'),
      easyClose = TRUE,size = "xl",
      footer = downloadButton("btn_dl_vcf_report_html","Download",icon = icon("download"),class="btn-primary")
    ))
  })
  
  output$btn_dl_vcf_report_html <- downloadHandler(
    filename = function(){str_glue("{basename(current_vcf())}.{input$resistance_db}.html")},
    content = function(file) {file.copy(str_glue("{current_vcf()}.{input$resistance_db}.html"),file)}
  )
  
  output$btn_dl_vcf_report_docx <- downloadHandler(
    filename = function(){str_glue("{basename(current_vcf())}.{input$resistance_db}.docx")},
    content = function(file) {file.copy(str_glue("{current_vcf()}.{input$resistance_db}.docx"),file)}
  )

}

shinyApp(ui = ui(),server = server)
