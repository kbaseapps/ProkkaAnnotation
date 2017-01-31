
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
 * <p>Original spec-file type: AnnotateContigsOutput</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "output_genome_ref"
})
public class AnnotateContigsOutput {

    @JsonProperty("output_genome_ref")
    private String outputGenomeRef;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("output_genome_ref")
    public String getOutputGenomeRef() {
        return outputGenomeRef;
    }

    @JsonProperty("output_genome_ref")
    public void setOutputGenomeRef(String outputGenomeRef) {
        this.outputGenomeRef = outputGenomeRef;
    }

    public AnnotateContigsOutput withOutputGenomeRef(String outputGenomeRef) {
        this.outputGenomeRef = outputGenomeRef;
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
        return ((((("AnnotateContigsOutput"+" [outputGenomeRef=")+ outputGenomeRef)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
