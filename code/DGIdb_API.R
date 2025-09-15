library(ghql)
library(jsonlite)
library(dplyr)

con <- GraphqlClient$new(
  url = "https://dgidb.org/api/graphql"
)  

gene_name <- "FLT1"

#####
query_text <- paste0('
{
  genes(names: ["',gene_name,'"]) {
    nodes {
      interactions {
        drug {
          name
          conceptId
        }
        interactionScore
        interactionTypes {
          type
          directionality
        }
        interactionAttributes {
          name
          value
        }
        publications {
          pmid
        }
        sources {
          sourceDbName
        }
      }
    }
  }
}
')

qry <- Query$new()
qry$query(name = "geneInteractions", x = query_text)

#####
res <- con$exec(qry$queries$geneInteractions)


#####
res_json <- fromJSON(res, flatten = TRUE)

#####
interactions <- res_json$data$genes$nodes
interactions1 <- interactions[[1]]
result <- interactions1[[1]]







