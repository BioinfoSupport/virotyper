library(bslib)
library(shiny)
library(bsicons)


options(shiny.maxRequestSize=5000*1024^2)


ui_vcf_panel <- function() {
  nav_panel(
    title = "VCF",
    tags$p("Use this form to generate a resistance report from a VCF file."),
    inputPanel(
      title = "Generate a report from a VCF containing variants",
      selectInput("resistance_db","Resistance database",choices = c("HSV1"="data/HHV1/db","HSV2"="data/HHV2/db")),
      fileInput("vcf_file","VCF file",accept = c(".vcf",".gz"),multiple = FALSE,placeholder = "Upload a VCF file (with INFO fields DP,AF,BCSQ)"),
    ),
    shinyjs::disabled(actionButton("btn_make_vcf_report","Generate VCF Report",icon = icon("cog"))),
    shinyjs::disabled(actionButton("btn_view_vcf_report","View HTML report",icon = icon("eye"))),
    shinyjs::disabled(downloadButton("btn_dl_vcf_report_docx","Download DOCX Report",icon = icon("download")))
  )
}

ui_bam_panel <- function() {
  nav_panel(
    title = "BAM",
    tags$p("Use this form to generate a report from a BAM file containing ONT reads aligned on a reference"),
    inputPanel(
      fileInput("bam_file","BAM file",accept = ".bam",multiple = FALSE,placeholder = "Upload BAM file of aligned reads on the reference",width="100%")
    ),
    shinyjs::disabled(actionButton("btn_make_bam_report","Generate BAM Report",icon = icon("cog"))),
    shinyjs::disabled(actionButton("btn_view_bam_report","View BAM Report",icon = icon("eye")))
  )
}

ui_fasta_panel <- function() {
  nav_panel(
    title="FASTA",
    markdown("Use this form to detect mutations in a consensus sequence relatively to a reference database.
                This app align the consensus sequences in the given FASTA file on the reference sequence using `minimap2`, and 
                detect the variants with `bcftools`."),
    inputPanel(
      tags$p(tags$b("Not available yet !!!"))
    )
  )
}

ui <- function() {
  page_navbar(
    title = icon("biohazard"),
    theme = bs_theme(),
    collapsible = TRUE,
    inverse = TRUE,
    header = shinyjs::useShinyjs(),
    footer = card(max_height="200px",
                  class="bg-warning",
                  card_header("Disclaimer"),
                  card_body("IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
          ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
          TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
          OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
                  )
    ),
    ui_vcf_panel(),
    ui_bam_panel(),
    ui_fasta_panel()
    #   selectInput("resistance_db","Resistance database",choices = c("HSV1","HSV2")),
    #   textAreaInput("fasta_seq","DNA sequences (FASTA)",placeholder = ">UL23\nATGCTCATGACTGCATCTA...\n>UL30\nATGCTCATGACTGCATCTA...",height = "400px"),
    #   shinycssloaders::withSpinner(textOutput("terminal_stdout",container = pre),type = 1)
  )
}

server <- function(input, output, session) {

  #-#-#-#-#-#-#-#
  # BAM
  #-#-#-#-#-#-#-#
  bam_report_dir <- reactiveVal(NULL)
  
  observeEvent(input$bam_file,{
    validate(need(is.character(input$bam_file$datapath),"invalid BAM"))
    validate(need(file.exists(input$bam_file$datapath),"invalid BAM"))
    shinyjs::enable("btn_make_bam_report")
  })
  
  observeEvent(input$btn_make_bam_report,{
    validate(need(is.character(input$bam_file$datapath),"invalid BAM"))
    shinyjs::disable("btn_view_bam_report")
    shinyjs::disable("btn_dl_bam_report_html")
    
    # Copy inputs to temp location
    tempdir <- tempdir()
    file.copy("notebooks/bam_report.Rmd", file.path(tempdir,"template.Rmd"), overwrite = TRUE)
    file.copy(input$bam_file$datapath, file.path(tempdir,input$bam_file$name), overwrite = TRUE)
    
    # Generate report
    withProgress(message="Generate BAM report",{
      rmarkdown::render(
        file.path(tempdir,"template.Rmd"),
        output_file = file.path(tempdir,"bam_report.html"),
        params = list(input_bam_file=file.path(tempdir,input$bam_file$name)),
        envir = new.env(parent = globalenv())
      )
    })
    bam_report_dir(tempdir)
    shinyjs::enable("btn_view_bam_report")
    shinyjs::enable("btn_dl_bam_report_html")
  })
  
  observeEvent(input$btn_view_bam_report,{
    showModal(modalDialog(
      title = "BAM report",
      tags$iframe(srcdoc=readLines(file.path(bam_report_dir(),"bam_report.html")), style='width:100%;height:600px;'),
      easyClose = TRUE,size = "xl",
      footer = downloadButton("btn_dl_bam_report_html","Download",icon = icon("download"),class="btn-primary")
    ))
  })
  
  output$btn_dl_bam_report_html <- downloadHandler(
    filename = function(){paste0(input$bam_file$name,".html")},
    content = function(file) {file.copy(file.path(bam_report_dir(),"bam_report.html"),file)}
  )
  
  #-#-#-#-#-#-#-#
  # VCF
  #-#-#-#-#-#-#-#
  vcf_report_dir <- reactiveVal(NULL)
  
  observeEvent(input$vcf_file,{
    validate(need(is.character(input$vcf_file$datapath),"invalid VCF"))
    shinyjs::enable("btn_make_vcf_report")
  })
  
  observeEvent(input$btn_make_vcf_report,{
    validate(need(is.character(input$vcf_file$datapath),"invalid VCF"))
    shinyjs::disable("btn_view_vcf_report")
    shinyjs::disable("btn_dl_vcf_report_docx")
    
    # Copy inputs to temp location
    tempdir <- tempdir()
    file.copy("notebooks/vcf_report.Rmd", file.path(tempdir,"template.Rmd"), overwrite = TRUE)
    file.copy(input$vcf_file$datapath, file.path(tempdir,input$vcf_file$name), overwrite = TRUE)
    file.copy(input$resistance_db, tempdir, overwrite = TRUE,recursive = TRUE)
    
    # Generate report
    withProgress(message="Generate VCF report",{
      rmarkdown::render(
        file.path(tempdir,"template.Rmd"),
        output_file = file.path(tempdir,"vcf_report.html"),
        params = list(
          input_vcf_file = file.path(tempdir,input$vcf_file$name),
          input_db_dir = file.path(tempdir,"db"),
          output_docx_report = file.path(tempdir,"vcf_report.docx")),
        envir = new.env(parent = globalenv())
      )
    })
    shinyjs::enable("btn_view_vcf_report")
    shinyjs::enable("btn_dl_vcf_report_docx")
    bam_report_dir(tempdir)
  })
  
  observeEvent(input$btn_view_vcf_report,{
    showModal(modalDialog(
      title = "VCF report",
      tags$iframe(srcdoc=readLines(file.path(bam_report_dir(),"vcf_report.html")), style='width:100%;height:600px;'),
      easyClose = TRUE,size = "xl",
      footer = downloadButton("btn_dl_vcf_report_html","Download",icon = icon("download"),class="btn-primary")
    ))
  })
  
  output$btn_dl_vcf_report_html <- downloadHandler(
    filename = function(){paste0(input$vcf_file$name,".html")},
    content = function(file) {file.copy(file.path(bam_report_dir(),"vcf_report.html"),file)}
  )
  
  output$btn_dl_vcf_report_docx <- downloadHandler(
    filename = function(){paste0(input$vcf_file$name,".docx")},
    content = function(file) {file.copy(file.path(bam_report_dir(),"vcf_report.docx"),file)}
  )
  
  # output$terminal_stdout <- eventReactive(input$make_doc,{
  #   validate(need(is.character(input$vcf_file$datapath),"invalid fasta"))
  #   cmd <- sprintf("make %s 2>&1",input$vcf_file$datapath)
  #   print(cmd)
  #   withProgress(message="Generate report",{
  #     pipe(cmd) |> readLines() |> paste(collapse = "\n")
  #   })
  # })
}

shinyApp(ui = ui(),server = server)
