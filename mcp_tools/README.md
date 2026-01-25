# MCP servers

To facilitate obtaining the most robust LLM-based evaluation results via the scoring prompt template in our manuscript, we developed mcp servers for deterministic tasks such as [**id mapping**](https://github.com/tf1993614/SuperSCC/blob/main/mcp_tools/id_mapping.py), [**intersection statistics calculation**](https://github.com/tf1993614/SuperSCC/blob/main/mcp_tools/calculate_intersection.py) and [**Fisher's exact test**](https://github.com/tf1993614/SuperSCC/blob/main/mcp_tools/fisher_test.py). 

To apply those mcp servers, you can follow the [**MCP official tutorial**](https://modelcontextprotocol.io/docs/develop/build-server#set-up-your-environment). 

For example, to register the id mapping mcp sever in Cline extension in VsCode, you should:

		# Create a new directory for id_mapping project
		uv init id_mapping
		cd id_mapping

		# Create virtual environment and activate it
		uv venv
		source .venv/bin/activate

		# Install dependencies
		uv add "mcp[cli]" httpx mygene

		# Then put the id_mapping.py file into id_mapping folder

Similaryly, for calculate_intersection and fisher_test mcp servers, you should:

		# Create a new directory for calculate_intersection project
		uv init calculate_intersection
		cd calculate_intersection

		# Create virtual environment and activate it
		uv venv
		source .venv/bin/activate

		# Install dependencies
		uv add "mcp[cli]" httpx 

		# Then put the calculate_intersection.py file into calculate_intersection folder


		# Create a new directory for fisher_test project
		uv init fisher_test
		cd fisher_test

		# Create virtual environment and activate it
		uv venv
		source .venv/bin/activate

		# Install dependencies
		uv add "mcp[cli]" httpx scipy

		# Then put the fisher_test.py file into fisher_test folder

Once above environment is set up, open mcp configuration panel in Cline extension in VsCode  

![img](https://github.com/tf1993614/SuperSCC/blob/main/mcp_tools/img/screenshot_1.png)


and then make a configure json file shown below:

		{
			"mcpServers": {
					"id_mapping": {
					"command": "uv",
					"args": [
						"--directory",
						"your path to id_mapping folder", # replace the place holder
						"run",
						"id_mapping.py"
					]
				},
					"fisher_test": {
					"command": "uv",
					"args": [
						"--directory",
						"your path to fisher_test folder", # replace the place holder
						"run",
						"fisher_test.py"
					]
				},
					"calculate_intersection": {
					"command": "uv",
					"args": [
						"--directory",
						"your path to calculate_intersection folder", # replace the place holder
						"run",
						"calculate_intersection.py"
					]
				}
			}
		}

Once it's done, your Cline will show those mcp servers like:

![img](https://github.com/tf1993614/SuperSCC/blob/main/mcp_tools/img/screenshot_2.png)


## Example
The example of our prompt serving as a high-level directive for the LLM to execute the necessary mcp servers shows below and full running record can be found [**here**](https://github.com/tf1993614/SuperSCC/blob/main/mcp_tools/cline_full_record.md):

![img](https://github.com/tf1993614/SuperSCC/blob/main/mcp_tools/img/screenshot_3.gif)

