```mermaid
sequenceDiagram

Frontend->>Storage: Save input files.
Storage->>Tool-Container: Script activates tool.
Tool-Container->>Storage: Tool finishes job, produces output files.
Storage->>Frontend: Generate download link / direct user to the output file.

```

**Frontend**: HTML/JS/CSS supported simple interface
**Storage**: Temporary/Permanant space to manage tool-based scripts, inputs, and outputs.
**Tool-Container**: Docker containers that will be hosted on HVP server, directly takes input and outputs from the storage directory.