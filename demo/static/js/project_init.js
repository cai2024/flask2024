document.getElementById('confirmBtn').addEventListener('click', confirmSteps);
var raw_myOptions = [localStorage.getItem('selectedFileType') || "raw.fq"];
var fq_type = localStorage.getItem('selectedFqType') || "single";
var myOptions = [];

$(document).ready(function() {
    $('input[name="fileType"][value="' + raw_myOptions[0] + '"]').prop('checked', true);
    myOptions = [raw_myOptions[0]];
    $('input[name="fileType"]').change(function() {
        raw_myOptions = [$(this).val()];
        myOptions = [raw_myOptions[0]];
        localStorage.setItem('selectedFileType', $(this).val());
        location.reload();
    });
    $('input[name="fqType"][value="' + fq_type + '"]').prop('checked', true);
    $('input[name="fqType"]').change(function() {
        fq_type = $(this).val();
        localStorage.setItem('selectedFqType', fq_type);
        location.reload();
    });
    

    
});







// 初始化
async function confirmSteps() {
    confirmBtn.disabled = true;  // 禁用按钮
    console.log('steps:', steps);
    for (let i = 0; i < steps.length; i++) {
        steps[i].id = i;
    }
    for (const step of steps) {
        await initOptions(step.type, step.id);
    }
    console.log('Before disable:', $('#selected-steps'));
    $('#selected-steps').sortable('disable');
    const submitBtn = document.getElementById('submitBtn');
    console.log('Before disable:', $('#selected-steps'));
    submitBtn.disabled = false; 
}


// 初始化当前行的input和options
async function initOptions(stepType, stepId) {
    try {
        const response = await fetch('/get_options', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ step_type: stepType, my_options: myOptions })
        });
        const data = await response.json();
        const options = data.options;  // 核心任务，传递给更新UI
        await initOptionsUI(stepId, stepType, options);

        // 检查 stepType 的前三个字符是否为 "map"
        if (stepType.slice(0, 3) === 'map') {
            const mapResponse = await fetch('/index_select', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ step_type: stepType })
            });
            const mapData = await mapResponse.json();
            const index_options = mapData.index_options;
            await initIndex(stepId, stepType, index_options);
        }
    } catch (error) {
        console.error('Error:', error);
    }
}


async function initOptionsUI(stepId, stepType, options) {
    const step = steps.find(step => step.id === stepId);
    const optionsHtml = options.map(option => `<option value="${option}">${option}</option>`).join('');
    const inputDiv = document.createElement('div');
    inputDiv.id = `input-step-${stepId}`;
    inputDiv.classList.add('form-group');
    inputDiv.innerHTML = `
        <label><strong>step${stepId}: ${stepType}:</strong></label>
        <input id="input-${stepId}"  type="text" class="form-control" placeholder="请输入内容" onchange="updateInput(${stepId})" value="${steps.find(step => step.id === stepId).params || ''}">
        <select id="select-${stepId}" class="form-control mt-2" onchange="updateSelectOption(${stepId})">${optionsHtml}</select>
    `;
    document.getElementById('inputArea').appendChild(inputDiv);
    // 绑定输入框
    const inputElement = document.querySelector(`#input-step-${stepId} input`);
    step.params = inputElement.value;
    inputElement.addEventListener('input', () => {
        updateOutputArea();
    });
    // 绑定选择框
    const selectElement = document.querySelector(`#input-step-${stepId} select`);
    if (selectElement && options.length > 0) {
        const longestIndex = options.reduce((maxIndex, option, index, array) => option.length > array[maxIndex].length ? index : maxIndex, 0);
        selectElement.selectedIndex = longestIndex;
        step.select_option = selectElement.value;
        await updateMyOptions(stepType, stepId);
    }
}

async function initIndex(stepId, stepType, index_options) {
    const step = steps.find(step => step.id === stepId);
    const optionsHtml = index_options.map(option => `<option value="${option}">${option}</option>`).join('');
    const inputDiv = document.createElement('div');
    inputDiv.id = `index-step-${stepId}`;
    inputDiv.classList.add('form-group');
    inputDiv.innerHTML = `
        <label><strong>step${stepId}: ${stepType} 参考基因组:</strong></label>
        <select class="form-control mt-2">${optionsHtml}</select>
    `;
    document.getElementById('indexArea').appendChild(inputDiv);
    const selectElement = inputDiv.querySelector('select');
    selectElement.onchange = () => {
        step.select_index=selectElement.value;
    };
    if (selectElement && index_options.length > 0) {
        step.select_index = selectElement.value;
    }
}



