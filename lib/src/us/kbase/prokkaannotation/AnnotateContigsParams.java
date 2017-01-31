
package us.kbase.prokkaannotation;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: AnnotateContigsParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "assembly_ref",
    "output_workspace",
    "output_genome_name"
})
public class AnnotateContigsParams {

    @JsonProperty("assembly_ref")
    private String assemblyRef;
    @JsonProperty("output_workspace")
    private String outputWorkspace;
    @JsonProperty("output_genome_name")
    private String outputGenomeName;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("assembly_ref")
    public String getAssemblyRef() {
        return assemblyRef;
    }

    @JsonProperty("assembly_ref")
    public void setAssemblyRef(String assemblyRef) {
        this.assemblyRef = assemblyRef;
    }

    public AnnotateContigsParams withAssemblyRef(String assemblyRef) {
        this.assemblyRef = assemblyRef;
        return this;
    }

    @JsonProperty("output_workspace")
    public String getOutputWorkspace() {
        return outputWorkspace;
    }

    @JsonProperty("output_workspace")
    public void setOutputWorkspace(String outputWorkspace) {
        this.outputWorkspace = outputWorkspace;
    }

    public AnnotateContigsParams withOutputWorkspace(String outputWorkspace) {
        this.outputWorkspace = outputWorkspace;
        return this;
    }

    @JsonProperty("output_genome_name")
    public String getOutputGenomeName() {
        return outputGenomeName;
    }

    @JsonProperty("output_genome_name")
    public void setOutputGenomeName(String outputGenomeName) {
        this.outputGenomeName = outputGenomeName;
    }

    public AnnotateContigsParams withOutputGenomeName(String outputGenomeName) {
        this.outputGenomeName = outputGenomeName;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((("AnnotateContigsParams"+" [assemblyRef=")+ assemblyRef)+", outputWorkspace=")+ outputWorkspace)+", outputGenomeName=")+ outputGenomeName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
