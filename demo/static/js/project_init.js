
document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
async function confirmSteps() {
    confirmBtn.disabled = true;  // 禁用确认按钮
    for (let i = 0; i < steps.length; i++) {steps[i].id = i;}
    try {
        const response = await fetch('/common/init_div_steps', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ steps :steps })});
        const data = await response.json();
        initOptionsUI(data);
        soft_step = data;
    } catch (error) { console.error('Error:', error); }
    
    $('#selected-steps').sortable('disable'); // 禁用左侧按钮卡排序
    document.getElementById('submitBtn').disabled = false; // 允许提交
}


async function initOptionsUI(data) {
    data.forEach((step, index) => {
        let html = `
            <div id="input-step-${step.id}" class="form-group">
                <label><strong>step${step.id}: ${step.type}:</strong></label>
                <input id="input-${step.id}" type="text" class="form-control" placeholder="输入额外参数（默认为空）" onchange="updateInput(${step.id})" oninput="updateOutputArea()" value="${step.params || ''}">
        `;
        Object.entries(step.selectableOptions).forEach(([key, options]) => {
            html += `<select id="select-${step.id}-${key}" class="form-control mt-2" onchange='updateSelectOption("${step.id}", "${key}")' title="${key}">`;
            options.forEach(option => {
                const isSelected = step.rec_option && step.rec_option[key] === option;
                html += `<option value="${option}" ${isSelected ? 'selected' : ''}>${option}</option>`;
                
               
            });
            html += `</select>`;
        });
        html += `</div>`;

        if (step.type.slice(0, 3) === 'map' || step.type === 'BAM CGS_call_modify' || step.type === 'BAM TBS_call_modify') {
            let index_html = `
                <div id="index-step-${step.id}" class="form-group"> 
                <label><strong>step${step.id}: ${step.type} 参考基因组:</strong></label>
                
            `
            Object.entries(step.selectableIndex).forEach(([key, options]) => {
                index_html += `<select id="index-${step.id}-${key}" class="form-control mt-2" title="${key}">`;
                options.forEach(option => {
                    index_html += `<option value="${option}">${option}</option>`;
                });
                index_html += `</select>`;
            });
            index_html += `</div>`;
            document.getElementById('indexArea').insertAdjacentHTML('beforeend', index_html);
        }

        document.getElementById('inputArea').insertAdjacentHTML('beforeend', html);
        

    })
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


