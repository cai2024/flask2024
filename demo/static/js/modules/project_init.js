

function escapeHTML(str) {
    if (!str) return '';
    return str.replace(/&/g, "&amp;")
              .replace(/</g, "&lt;")
              .replace(/>/g, "&gt;")
              .replace(/"/g, "&quot;")
              .replace(/'/g, "&#039;");
}


async function initOptionsUI(data) {
    data.forEach((step, index) => {
        let html = `
            <div id="input-step-${step.id}" class="form-group">
                <label><strong>step${step.id}: ${step.type}:</strong></label>
                <input id="input-${step.id}" type="text" class="form-control" placeholder="输入额外参数（默认为空）"  value="${escapeHTML(step.params || '')}">
        `;
        Object.entries(step.selectableOptions).forEach(([key, options]) => {
            html += `<select id="select-${step.id}-${key}" class="form-control mt-2" onchange='updateSelectOption("${step.id}", "${key}")' title="${key}">`;

            options.forEach((option, index) => {
                const isSelected = step.rec_option[key] >= 0 && step.rec_option[key] < options.length && step.rec_option[key] === index;
                html += `<option value="${option}" ${isSelected ? 'selected' : ''}>${option}</option>`;
            });
            html += `</select>`;
        });
        html += `</div>`;

        if (step.type.slice(0, 3) === 'map' || step.type === 'BAM CGS_call_modify' || step.type === 'BAM TBS_call_modify' || step.type === 'BAM CT_call_modify' || step.type === 'BAM CT_XM_bam' ||  step.type === 'BAM astair_call_modify') {
            let index_html = `<div id="index-step-${step.id}" class="form-group"> <label><strong>step${step.id}: ${step.type} 参考基因组:</strong></label>`
            Object.entries(step.selectableIndex).forEach(([key, options]) => {
                if (key === "ref") {
                    index_html += `<select id="species-${step.id}" class="form-control mt-2" title="step${step.id}_ref_species" onchange="handleSpeciesChange(this, ${step.id})">`;
                    Object.keys(options).forEach(species => {
                    index_html += `<option value="${species}">${species.charAt(0).toUpperCase() + species.slice(1)}</option>`;});
                    index_html += `</select> <select id="version-${step.id}" class="form-control mt-2" title="step${step.id}_ref_release"> </select>`;
                }else {
                    index_html += `<select id="index-${step.id}-${key}" class="form-control mt-2" title="${key}">`;
                    options.forEach(option => {index_html += `<option value="${option}">${option}</option>`;});
                    index_html += `</select>`;
                }
            });
            index_html += `</div>`;
            document.getElementById('indexArea').insertAdjacentHTML('beforeend', index_html);
            handleSpeciesChange(document.getElementById(`species-${step.id}`) ,step.id);
        }
        document.getElementById('inputArea').insertAdjacentHTML('beforeend', html);
    })
}

function handleSpeciesChange(selectElement,stepId) {
    const step = steps[stepId]
    const selectedSpecies = selectElement.value;
    const versionSelect = document.getElementById(`version-${step.id}`);
    versionSelect.innerHTML = '';
    if (selectedSpecies && step.selectableIndex["ref"][selectedSpecies]) {
        step.selectableIndex["ref"][selectedSpecies].forEach(version => {
            const option = document.createElement("option");
            option.value = version;
            option.textContent = version;
            versionSelect.appendChild(option);
        });
    }
}



async  function updateSelectOption(stepId, change_key) {
    const step = steps.find(step => step.id === stepId);
    try {
        const selectedValue = document.querySelector(`#input-step-${stepId} select#select-${stepId}-${change_key}`).value
        if (change_key === "fq1") {
            const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-fq2`)
            if (selectElement) {
                const index = Array.from(selectElement.options).findIndex(option => option.value === selectedValue.replace(/_1\.fq\.gz$/, '_2.fq.gz'));
                if (index !== -1) {
                    selectElement.selectedIndex = index;
                }
            }
        } else if (change_key === "fq2") {
            const selectElement = document.querySelector(`#input-step-${stepId} select#select-${stepId}-fq1`)
            if (selectElement) {
                const index = Array.from(selectElement.options).findIndex(option => option.value === selectedValue.replace(/_2\.fq\.gz$/, '_1.fq.gz'));
                if (index !== -1) {
                    selectElement.selectedIndex = index;
                }
            }
        }
    }catch (error) { console.error('Error:', error); }
}



function collectStepInputs(steps) {
    steps.forEach((step, stepId) => {
        step.params = document.getElementById(`input-${stepId}`).value;
        const selectElements = document.querySelectorAll(`#input-step-${stepId} select`);
        selectElements.forEach(select => {step.select_option[select.title] = select.value;});
        const indexElements = document.querySelectorAll(`#index-step-${stepId} select`);
        indexElements.forEach(select => {step.select_index[select.title] = select.value;});
    });
}

