library(bs4Dash)


ui <- dashboardPage(
  dark = NULL,help = NULL,
  header = dashboardHeader(icon("biohazard")," Viro Resistance Typer"),
  sidebar = dashboardSidebar(disable=TRUE),
  body = dashboardBody(
    fluidRow(
      box(title="Disclaimer",width=12,collapsible = FALSE,
          "IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
          DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
          ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
          IN THE SOFTWARE.
      "),
      box(title="Input",width=12,collapsible = FALSE,
          selectInput("resistance_db","Resistance database",choices = c("HSV1","HSV2")),
          fileInput("vcf_file","VCF file",placeholder = "Upload a VCF file (with INFO fields DP,AF,BCSQ)"),
          fileInput("bam_file","BAM file",placeholder = "Upload corresponding BAM file (optional)"),
          #textAreaInput("fasta_seq","DNA sequences (FASTA)",placeholder = ">UL23\nATGCTCATGACTGCATCTA...\n>UL30\nATGCTCATGACTGCATCTA...",height = "400px"),
          actionButton("make_doc","Report.doc",icon = icon("download")),
          shinycssloaders::withSpinner(textOutput("terminal_stdout",container = pre),type = 1)
      )
    )
  )
)

server <- function(input, output, session) {
  output$terminal_stdout <- eventReactive(input$make_doc,{
    validate(need(is.character(input$vcf_file$datapath),"invalid fasta"))
    cmd <- sprintf("make %s 2>&1",input$vcf_file$datapath)
    print(cmd)
    withProgress(message="Generate report",{
      pipe(cmd) |> readLines() |> paste(collapse = "\n")
    })
  })
}

shinyApp(ui, server)
