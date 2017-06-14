import aiohttp
import asyncio

# 
# ### Pull data from cytobank.

# cytobank_request <- function (..., req) {
#     cy_sesh <- authenticate(...)
#     tryCatch({
#         match.fun(req)(cy_sesh)
#     },
#     finally = authentication.logout(cy_sesh))
# }

# download_all_fcs <- function (auth, exp_id, dir = getwd(),
#                               n = pretty_good_parallelism()) {
#     fcs_listings <- fcs_files.list(auth, experiment_id = exp_id) %>%
#         select(id, filename, md5sum)
#     ## throw unless we can download everything
#     ## return downloaded file paths
# }

# dl_check_fcs <- function (auth, exp_id, fcs_id, file_md5) {
#     fcs_dl <- fcs_files.download(
#         auth, experiment_id = exp_id, fcs_file_id = fcs_id)

# }

# get_gates_pops_set <- function (auth, exp_id, dir = getwd(),
#                                 without_fcs_dl = FALSE) {
#     ## download every fcs file from the experiment
#     if (!without_fcs_dl) {
#         all_fcs_zipped <- fcs_files.download_zip(auth, experiment_id = exp_id)
#         unzip(all)
#     }
#     ## download gatingml as xml
#     gating_file <- gates.gatingML_download(auth, experiment_id = exp_id)
#     cytobank2GatingSet(gating_file, list.files(path = "./5-2"))
# }

import aiohttp
import asyncio

async def stream_to_file(stream, path, mode = 'wb'):
    while

async def fetch(client):
    async with client.get('http://python.org') as resp:
        assert resp.status == 200
        return await resp.text()

async def main(loop):
    async with aiohttp.ClientSession(loop=loop) as client:
        html = await fetch(client)
        print(html)

loop = asyncio.get_event_loop()
loop.run_until_complete(main(loop))
